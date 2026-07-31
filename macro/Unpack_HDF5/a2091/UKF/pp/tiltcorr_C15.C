/// @file tiltcorr_C15.C
/// @brief Does straightening the Ex-vs-angle tilt improve the Ex resolution, and is it safe?
///
/// A tilted locus projected onto the Ex axis IS wider than it needs to be, so straightening it is
/// a legitimate way to recover resolution. The question is which variable to straighten in.
///
///   theta_lab / vertex z : measured INDEPENDENTLY of Ex. An instrumental bias (path length,
///        drift, angle) is a function of these, so subtracting a fitted a(theta_lab) removes real
///        smearing. Safe.
///   theta_cm            : computed FROM Ex inside kine_2b. Subtracting g(theta_cm) subtracts a
///        function of each event's OWN Ex, which can compress the distribution self-referentially:
///        the peak narrows while LEVEL SPACINGS distort. That is not resolution, it is collapse.
///
/// 15C(p,d)14C carries its own ruler -- the 14C 6.09/6.59/6.73/7.01 cluster -- so the two schemes
/// can be judged on more than the g.s. width:
///   * a correction that improves resolution keeps the g.s.-to-cluster SPACING intact
///   * a correction that compresses shrinks that spacing while narrowing everything
///
/// The g.s. is isolated on its KINEMATIC LOCUS, not with a flat Ex window: at each theta_lab the
/// g.s. sits at the KE solving Ex(KE) = 0. A flat window mixed in the 6-7 MeV cluster, which made
/// <KE> non-monotonic in theta and produced absurd implied KE biases of 9-11 MeV.
///
///   root -b -q 'pp/tiltcorr_C15.C()'

static double om(double x, double y, double z) { return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

struct Kin { double ex, thcm; };
static Kin kine(double m1,double m2,double m3,double m4,double K,double th,double Ke)
{
   double E1=K+m1, E3=Ke+m3, E4=E1+m2-E3;
   double s=m1*m1+m2*m2+2*m2*E1, u=m2*m2+m3*m3-2*m2*E3;
   double a=(cos(th)*om(s,m1*m1,m2*m2)*om(u,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)+s+u-m2*m2;
   if(a<0) return {std::nan(""),std::nan("")};
   double m4x=sqrt(a), ex=m4x-m4, t=m2*m2+m4x*m4x-2*m2*E4;
   double g=(s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4x*m4x)+(m1*m1-m2*m2)*(m3*m3-m4x*m4x))
            /(om(s,m1*m1,m2*m2)*om(s,m3*m3,m4x*m4x));
   if(g<-1)g=-1; if(g>1)g=1;
   return {ex,(TMath::Pi()-acos(g))*TMath::RadToDeg()};
}

/// KE at which Ex(KE) = exTarget for this theta_lab (Ex falls monotonically with KE -> bisect)
static double keAtEx(double m1,double m2,double m3,double m4,double E,double th,double exTarget)
{
   double lo=0.2, hi=120.0;
   double flo=kine(m1,m2,m3,m4,E,th,lo).ex, fhi=kine(m1,m2,m3,m4,E,th,hi).ex;
   if(std::isnan(flo)||std::isnan(fhi)||(flo-exTarget)*(fhi-exTarget)>0) return std::nan("");
   for(int i=0;i<80;i++){
      double mid=0.5*(lo+hi), fm=kine(m1,m2,m3,m4,E,th,mid).ex;
      if(std::isnan(fm)) { hi=mid; continue; }
      if((flo-exTarget)*(fm-exTarget)<=0){ hi=mid; fhi=fm; } else { lo=mid; flo=fm; }
   }
   return 0.5*(lo+hi);
}

static bool fitPk(TH1F*h,double c,double win,double&mu,double&er,double&sg)
{
   if(h->GetEntries()<80) return false;
   TF1 g("g","gaus(0)+pol1(3)",c-win,c+win);
   int bm=h->GetMaximumBin();
   g.SetParameters(h->GetBinContent(bm),c,0.35,0,0);
   g.SetParLimits(1,c-win,c+win); g.SetParLimits(2,0.03,1.5);
   if(h->Fit(&g,"QRN")!=0) return false;
   mu=g.GetParameter(1); er=g.GetParError(1); sg=g.GetParameter(2);
   return er>0 && er<0.5;
}

void tiltcorr_C15(TString cache="plots/proton_kin_pd_ukf.root", double Ebeam=170.0,
                  double mEjectAmu=2.014102, double mResidAmu=14.003242,
                  double keWin=3.0, double clusterEx=6.5)
{
   gStyle->SetOptStat(0);
   const double u=931.49401, m1=15.0105993*u, m2=1.007825*u, m3=mEjectAmu*u, m4=mResidAmu*u;
   TString dir=gSystem->DirName(__FILE__);
   TFile*f=TFile::Open(cache.BeginsWith("/")?cache:dir+"/"+cache);
   TNtuple*n=(TNtuple*)f->Get("pk");
   float ke,th,vz,tc,ex,c2;
   n->SetBranchAddress("ke",&ke);n->SetBranchAddress("theta",&th);n->SetBranchAddress("vertexz",&vz);
   n->SetBranchAddress("thcm",&tc);n->SetBranchAddress("ex",&ex);n->SetBranchAddress("chi2ndf",&c2);

   struct Ev { double ke,th,vz,ex,thcm; };
   std::vector<Ev> ev;
   for(Long64_t i=0;i<n->GetEntries();i++){
      n->GetEntry(i); if(c2>5||ke<=0) continue;
      Kin k=kine(m1,m2,m3,m4,Ebeam,th*TMath::DegToRad(),ke);
      if(std::isnan(k.ex)) continue;
      ev.push_back({ke,th,vz,k.ex,k.thcm});
   }
   printf("\n=== tilt correction test: %s ===\ntracks %zu, Ebeam %.0f\n",
          gSystem->BaseName(cache),ev.size(),Ebeam);

   // ---- 1. g.s. offset vs theta_lab, isolated ON THE LOCUS ----
   printf("\n-- g.s. on its kinematic locus (|KE - KE_gs(theta)| < %.1f MeV) --\n",keWin);
   printf("%-12s %8s %9s %9s %8s %8s\n","theta_lab","N","KE_gs","mu_Ex","err","sigma");
   const int NB=10; double t0=10, dt=4;
   std::vector<double> bt, bmu;
   for(int b=0;b<NB;b++){
      double tlo=t0+b*dt, thi=tlo+dt, tc0=(tlo+thi)/2*TMath::DegToRad();
      double kegs=keAtEx(m1,m2,m3,m4,Ebeam,tc0,0.0);
      if(std::isnan(kegs)) continue;
      TH1F h("h","",300,-6,10);
      for(auto&e:ev) if(e.th>=tlo&&e.th<thi&&fabs(e.ke-kegs)<keWin) h.Fill(e.ex);
      double mu,er,sg;
      if(fitPk(&h,0.0,1.6,mu,er,sg)){
         printf("%-12s %8.0f %9.2f %9.3f %9.3f %8.3f\n",Form("%.0f-%.0f",tlo,thi),h.GetEntries(),kegs,mu,er,sg);
         bt.push_back((tlo+thi)/2); bmu.push_back(mu);
      } else printf("%-12s %8.0f %9.2f %9s\n",Form("%.0f-%.0f",tlo,thi),h.GetEntries(),kegs,"(no fit)");
   }
   if(bt.size()<3){ printf("too few locus slices\n"); return; }

   // linear a(theta_lab) through those offsets
   TGraph gl(bt.size()); for(size_t i=0;i<bt.size();i++) gl.SetPoint(i,bt[i],bmu[i]);
   TF1 la("la","pol1",bt.front()-2,bt.back()+2); gl.Fit(&la,"QRN");
   printf("  a(theta_lab) = %+.4f %+.4f*theta   (MeV)\n",la.GetParameter(0),la.GetParameter(1));

   // ---- 2. g.s. offset vs theta_cm, same locus selection ----
   std::vector<double> ct, cmu;
   for(double c=20;c<110;c+=12){
      TH1F h("h","",300,-6,10);
      for(auto&e:ev){
         double tr=e.th*TMath::DegToRad();
         double kegs=keAtEx(m1,m2,m3,m4,Ebeam,tr,0.0);
         if(std::isnan(kegs)||fabs(e.ke-kegs)>keWin) continue;
         if(e.thcm>=c&&e.thcm<c+12) h.Fill(e.ex);
      }
      double mu,er,sg;
      if(fitPk(&h,0.0,1.6,mu,er,sg)){ ct.push_back(c+6); cmu.push_back(mu); }
   }
   TF1 lc("lc","pol1",0,180);
   if(ct.size()>=3){ TGraph gc(ct.size());
      for(size_t i=0;i<ct.size();i++) gc.SetPoint(i,ct[i],cmu[i]);
      gc.Fit(&lc,"QRN");
      printf("  g(theta_cm)  = %+.4f %+.4f*theta_cm (MeV)\n",lc.GetParameter(0),lc.GetParameter(1)); }

   // ---- 3. apply each correction and measure g.s. AND the 14C cluster ----
   printf("\n-- effect of each correction (all events, no locus cut) --\n");
   printf("%-22s %9s %9s %9s %9s %11s\n","scheme","gs_mu","gs_sigma","clu_mu","clu_sigma","spacing");
   auto report=[&](const char*name,int mode){
      TH1F h("h","",400,-8,16);
      for(auto&e:ev){
         double c=0;
         if(mode==1) c=la.Eval(e.th);
         if(mode==2) c=lc.Eval(e.thcm);
         h.Fill(e.ex-c);
      }
      double g1,g1e,g1s, c1,c1e,c1s;
      bool og=fitPk(&h,0.0,1.6,g1,g1e,g1s);
      bool oc=fitPk(&h,clusterEx,1.8,c1,c1e,c1s);
      printf("%-22s %9s %9s %9s %9s %11s\n",name,
             og?Form("%+.3f",g1):"-", og?Form("%.3f",g1s):"-",
             oc?Form("%+.3f",c1):"-", oc?Form("%.3f",c1s):"-",
             (og&&oc)?Form("%.3f",c1-g1):"-");
      return h.GetEntries();
   };
   report("none",0);
   report("minus a(theta_lab)",1);
   report("minus g(theta_cm)",2);
   printf("\nJudge on SPACING, not width: the true 14C cluster sits ~6.1-7.0 MeV above the g.s.\n");
   printf("A scheme that narrows the peaks but shrinks the spacing is compressing the spectrum.\n");

   TCanvas*c=new TCanvas("ctc","tilt corr",1500,520); c->Divide(3,1);
   const char*nm[3]={"none","minus a(#theta_{lab})","minus g(#theta_{cm})"};
   for(int mode=0;mode<3;mode++){
      TH1F*h=new TH1F(Form("hh%d",mode),Form("%s;E_{x}(^{14}C) [MeV];counts",nm[mode]),300,-8,16);
      for(auto&e:ev){ double cc=0; if(mode==1)cc=la.Eval(e.th); if(mode==2)cc=lc.Eval(e.thcm); h->Fill(e.ex-cc); }
      c->cd(mode+1); gPad->SetLogy(); h->Draw();
   }
   c->SaveAs(dir+"/plots/tiltcorr_pd.png");
   f->Close();
}
