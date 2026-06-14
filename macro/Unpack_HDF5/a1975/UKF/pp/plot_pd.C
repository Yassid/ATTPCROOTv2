/// @file plot_pd.C
/// @brief Static (PNG) 4-panel visualization of the 16C(p,d)15C deuteron production —
/// the headless counterpart of the interactive explore_pd.C GUI. Reads the cached
/// kinematics ntuple (/tmp/pd_kin.root from pp/cache_pd_run.C) and draws:
///   (1) 15C E_x spectrum with a Gaussian g.s. fit + 15C state markers,
///   (2) deuteron KE vs theta_lab with 2-body kinematic lines,
///   (3) E_x vs theta_cm,  (4) vertex z vs E_x.
/// Works over X11/headless (no Eve/OpenGL); writes a PNG.
///
///   root -l -b -q 'pp/plot_pd.C'                       // defaults: /tmp/pd_kin.root, Ebeam=192
///   root -l -b -q 'pp/plot_pd.C("/tmp/pd_kin.root",192,5,"/tmp/pd_spectrum.png")'

static double pd_om2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
// returns {Ex(15C), theta_cm[deg]} for p(16C,15C*)d (relativistic 2-body)
static std::pair<double,double> pd_kine(double Kp,double thl,double Ke){
   const double u=931.49401; double m1=16.0147*u,m2=1.007825*u,m3=2.014102*u,m4=15.010599*u;
   double Et1=Kp+m1,Et3=Ke+m3,Et4=Et1+m2-Et3; double s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double a=(cos(thl)*pd_om2(s,m1*m1,m2*m2)*pd_om2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   if(a<0) return {NAN,NAN}; double m4x=std::sqrt(a); double ex=m4x-m4;
   double t=m2*m2+m4x*m4x-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4x*m4x)+(m1*m1-m2*m2)*(m3*m3-m4x*m4x))
                                    /(pd_om2(s,m1*m1,m2*m2)*pd_om2(s,m3*m3,m4x*m4x)));
   return {ex, tcm*TMath::RadToDeg()};
}
// deuteron KE-vs-theta_lab locus for 15C excitation Ex
static TGraph *pd_line(double Eb,double Ex,Color_t col,int style){
   const double u=931.49401; double m1=16.0147*u,m2=1.007825*u,m3=2.014102*u,m4=15.010599*u+Ex;
   double E1=Eb+m1,P=std::sqrt(E1*E1-m1*m1),W=E1+m2,s=W*W-P*P;
   auto*g=new TGraph(); g->SetLineColor(col); g->SetLineWidth(2); g->SetLineStyle(style);
   if(s<=(m3+m4)*(m3+m4)) return g;
   double rs=std::sqrt(s),E3cm=(s+m3*m3-m4*m4)/(2*rs),p3cm=std::sqrt(std::max(0.,E3cm*E3cm-m3*m3));
   double beta=P/W,gamma=W/rs;
   for(double tc=0;tc<=180;tc+=0.5){ double c=std::cos(tc*TMath::DegToRad()),sn=std::sin(tc*TMath::DegToRad());
      double Elab=gamma*(E3cm+beta*p3cm*c),pz=gamma*(p3cm*c+beta*E3cm),pp=p3cm*sn;
      double th=std::atan2(pp,pz)*TMath::RadToDeg(),ke=Elab-m3; if(ke>0&&th>=0) g->SetPoint(g->GetN(),th,ke); }
   return g;
}

void plot_pd(TString cache="/tmp/pd_kin.root", double Ebeam=192, double chi2max=5,
             double icMin=950, double icMax=1350, TString out="/tmp/pd_spectrum.png")
{
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TFile*f=TFile::Open(cache);
   if(!f||f->IsZombie()){ printf("cannot open %s (run pp/cache_pd_run.C first)\n",cache.Data()); return; }
   TTree*t=(TTree*)f->Get("pk"); if(!t){ printf("no tree pk\n"); return; }
   float ke,theta,vz,chi2ndf,ic; t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta);
   t->SetBranchAddress("vz",&vz); t->SetBranchAddress("chi2ndf",&chi2ndf); t->SetBranchAddress("ic",&ic);

   auto*hEx=new TH1F("hEx",Form("^{16}C(p,d)^{15}C  E_{x}(^{15}C), E_{beam}=%.0f;E_{x} [MeV];deuterons",Ebeam),160,-5,12);
   auto*hKT=new TH2F("hKT","deuteron KE vs #theta_{lab};#theta_{lab} [deg];KE [MeV]",120,0,95,120,0,45);
   auto*hEt=new TH2F("hEt","E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]",90,0,180,160,-5,12);
   auto*hZE=new TH2F("hZE","vertex z vs E_{x};E_{x} [MeV];vertex z [mm]",160,-5,12,110,-50,1000);

   long n=0;
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(chi2ndf>chi2max) continue;
      if(ic<icMin||ic>icMax) continue; // ion-chamber 16C beam gate (set icMin<0 to disable)
      hKT->Fill(theta,ke); auto[ex,tcm]=pd_kine(Ebeam,theta*TMath::DegToRad(),ke);
      if(!std::isnan(ex)){ hEx->Fill(ex); hEt->Fill(tcm,ex); hZE->Fill(ex,vz); } ++n; }

   double mn=0,sg=0;
   if(hEx->GetEntries()>50){ hEx->Fit("gaus","Q0","",-2,2);
      if(hEx->GetFunction("gaus")){ mn=hEx->GetFunction("gaus")->GetParameter(1); sg=hEx->GetFunction("gaus")->GetParameter(2); } }

   auto*c=new TCanvas("c_pd","16C(p,d)15C",1500,900); c->Divide(2,2);
   // 1) Ex spectrum + g.s. fit + 15C level markers
   c->cd(1); hEx->SetLineColor(kBlue+1); hEx->Draw("hist");
   if(hEx->GetFunction("gaus")){ hEx->GetFunction("gaus")->SetNpx(500); hEx->GetFunction("gaus")->SetLineColor(kRed+1); hEx->GetFunction("gaus")->Draw("same"); }
   { double ym=hEx->GetMaximum(); struct S{double e;const char*n;}; S st[4]={{0,"g.s."},{0.74,"0.74"},{3.10,"3.10"},{4.66,"4.66"}};
     auto*tx=new TLatex(); tx->SetTextSize(0.028); tx->SetTextAngle(90);
     for(auto&s:st){ auto*l=new TLine(s.e,0,s.e,ym); l->SetLineColor(s.e==0?kRed+1:kGray+2); l->SetLineStyle(2); l->Draw();
        tx->SetTextColor(s.e==0?kRed+1:kGray+2); tx->DrawLatex(s.e,0.4*ym,s.n); }
     auto*t2=new TLatex(); t2->SetNDC(); t2->SetTextSize(0.04); t2->DrawLatex(0.5,0.82,Form("g.s.: %.2f#pm%.2f MeV",mn,sg));
     t2->DrawLatex(0.5,0.76,Form("N=%ld deuterons",n)); }
   // 2) KE vs theta_lab + kinematic lines
   c->cd(2); hKT->Draw("colz");
   for(double ex:{0.0,3.10,4.66}){ auto*g=pd_line(Ebeam,ex,ex==0?kRed+1:kGray+2,ex==0?1:2); if(g->GetN()) g->Draw("L same"); }
   // 3) Ex vs theta_cm  4) vertex z vs Ex
   c->cd(3); hEt->Draw("colz");
   c->cd(4); hZE->Draw("colz");
   c->SaveAs(out);
   printf("16C(p,d)15C: %ld deuterons | g.s. peak %.3f +/- %.3f MeV | wrote %s\n",n,mn,sg,out.Data());
}
