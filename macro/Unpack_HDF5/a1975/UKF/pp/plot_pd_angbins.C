/// @file plot_pd_angbins.C
/// @brief 16C(p,d)15C excitation-energy spectra in theta_cm slices (default 10-deg
/// wide, starting at 10 deg) — one Ex panel per CM-angle bin, for angular-distribution
/// / DWBA work. Reads the cached deuteron kinematics (/tmp/pd_kin.root from
/// pp/cache_pd_run.C), applies the IC beam gate + chi2 cut, and draws a grid of Ex
/// histograms with the 15C level markers. Headless (PNG).
///
///   root -l -b -q 'pp/plot_pd_angbins.C'                                  // 10-deg bins from 10, 8 bins
///   root -l -b -q 'pp/plot_pd_angbins.C(10,10,9,-0.38)'                   // start,width,nbins,exShift

static double pd_om2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
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

void plot_pd_angbins(double thcmStart=10, double thcmWidth=10, int nbins=8, double exShift=-0.38,
                     double Ebeam=192, double chi2max=5, double icMin=950, double icMax=1350,
                     double exLo=-3, double exHi=11, int nex=140, TString cache="/tmp/pd_kin.root",
                     TString out="/tmp/pd_angbins.png")
{
   gStyle->SetOptStat(0);
   TFile*f=TFile::Open(cache);
   if(!f||f->IsZombie()){ printf("cannot open %s (run pp/cache_pd_run.C first)\n",cache.Data()); return; }
   TTree*t=(TTree*)f->Get("pk"); if(!t){ printf("no tree pk\n"); return; }
   float ke,theta,vz,chi2ndf,ic; t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta);
   t->SetBranchAddress("vz",&vz); t->SetBranchAddress("chi2ndf",&chi2ndf); t->SetBranchAddress("ic",&ic);

   std::vector<TH1F*> h(nbins);
   for(int k=0;k<nbins;++k){ double lo=thcmStart+k*thcmWidth, hi=lo+thcmWidth;
      h[k]=new TH1F(Form("hex%d",k),Form("#theta_{cm} = %.0f-%.0f#circ;E_{x}(^{15}C) [MeV];deuterons",lo,hi),nex,exLo,exHi); }

   std::vector<long> cnt(nbins,0);
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(chi2ndf>chi2max) continue; if(ic<icMin||ic>icMax) continue;
      auto[ex,tcm]=pd_kine(Ebeam,theta*TMath::DegToRad(),ke); if(std::isnan(ex)) continue;
      int k=(int)((tcm-thcmStart)/thcmWidth); if(k<0||k>=nbins) continue;
      h[k]->Fill(ex+exShift); cnt[k]++; }

   int nc=(nbins<=4)?nbins:((nbins<=6)?3:4), nr=(nbins+nc-1)/nc;
   auto*c=new TCanvas("c_ang","(p,d) Ex per theta_cm",380*nc,300*nr); c->Divide(nc,nr);
   const double ev[4]={0,0.74,3.10,4.66};
   for(int k=0;k<nbins;++k){ c->cd(k+1); h[k]->SetLineColor(kBlue+1); h[k]->Draw("hist");
      double ym=h[k]->GetMaximum(); if(ym<=0) continue;
      for(double e:ev){ auto*l=new TLine(e,0,e,ym); l->SetLineColor(e==0?kRed+1:kGray+2); l->SetLineStyle(2); l->Draw(); }
      double mn=0,sg=0; if(h[k]->GetEntries()>40){ h[k]->Fit("gaus","Q0","",-1.2,1.2);
         if(h[k]->GetFunction("gaus")){ mn=h[k]->GetFunction("gaus")->GetParameter(1); sg=h[k]->GetFunction("gaus")->GetParameter(2);
            h[k]->GetFunction("gaus")->SetLineColor(kRed+1); h[k]->GetFunction("gaus")->Draw("same"); } }
      auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.05);
      tx->DrawLatex(0.5,0.84,Form("N=%ld",cnt[k])); if(sg>0) tx->DrawLatex(0.5,0.77,Form("g.s.#sigma=%.2f",sg)); }
   c->SaveAs(out);
   printf("(p,d) Ex in %d theta_cm bins of %.0f deg from %.0f:",nbins,thcmWidth,thcmStart);
   for(int k=0;k<nbins;++k) printf(" [%.0f-%.0f]=%ld",thcmStart+k*thcmWidth,thcmStart+(k+1)*thcmWidth,cnt[k]);
   printf("\nwrote %s\n",out.Data());
}
