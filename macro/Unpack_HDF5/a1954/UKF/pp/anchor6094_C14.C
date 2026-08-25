/// @file anchor6094_C14.C
/// @brief Re-anchor the 14C excitation scale on the 6.094 MeV state for a given cache, then read
///        off where the 8.317 state lands and turn the difference into a residual gain on
///        (Ex - 6.094). This is the procedure of 12_excited.tex sec:exc-calib, re-run on the
///        2026-08-25 CATIMA production because the old anchor (Ebeam 164.25) belongs to the old
///        fitter and the refit moved the reconstructed energy by ~1.2 %.
///
///   root -b -q 'anchor6094_C14.C("plots/proton_kin_cmpB_cat.root","CATIMA matFX-on")'
static const double U = 931.49401;
static double om2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
static double ExOf(double Eb,double K,double thDeg)
{
   const double m1=14.003242*U,m2=1.007825*U,m3=1.007825*U,m4=14.003242*U;
   double Et1=Eb+m1,Et3=K+m3,s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double th=thDeg*TMath::DegToRad();
   double a=(std::cos(th)*om2(s,m1*m1,m2*m2)*om2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   return a<0 ? -99.0 : std::sqrt(a)-m4;
}
/// gaussian on a linear background over [lo,hi]; returns {mu, dmu, sigma, area, ok}
static std::tuple<double,double,double,double,bool>
FitPeak(TH1D *h,double lo,double hi,double sigSeed)
{
   int b0=h->FindBin(lo),b1=h->FindBin(hi);
   double best=-1; int bb=b0;
   for(int b=b0;b<=b1;++b) if(h->GetBinContent(b)>best){best=h->GetBinContent(b);bb=b;}
   double mu0=h->GetBinCenter(bb);
   TF1 f("f","gaus(0)+pol1(3)",lo,hi);
   f.SetParameters(best,mu0,sigSeed,0,0);
   f.SetParLimits(1,lo,hi);
   f.SetParLimits(2,0.05,1.2);          // bounds are REPORTED below, never silently trusted
   TFitResultPtr r=h->Fit(&f,"RQNS");
   bool ok = r.Get() && r->IsValid();
   double mu=f.GetParameter(1),sg=f.GetParameter(2);
   // rail check: a parameter sitting on its bound is a failure wearing a number
   if(std::fabs(sg-0.05)<1e-3||std::fabs(sg-1.2)<1e-3||std::fabs(mu-lo)<1e-3||std::fabs(mu-hi)<1e-3) ok=false;
   return {mu,f.GetParError(1),sg,f.GetParameter(0)*sg*std::sqrt(2*TMath::Pi())/h->GetBinWidth(1),ok};
}
void anchor6094_C14(TString cache="plots/proton_kin_cmpB_cat.root", TString label="CATIMA matFX-on",
                    double chi2Cut=5.0, double vzLo=20, double vzHi=1100,
                    double kcDenom=6428.6, double kcPivot=0.0, bool kcOn=true)
{
   gStyle->SetOptStat(0);
   TFile *f=TFile::Open(cache); if(!f||f->IsZombie()){printf("no %s\n",cache.Data());return;}
   TNtuple *t=(TNtuple*)f->Get("pk"); if(!t){printf("no pk tree\n");return;}
   std::vector<double> KE,TH_;
   float *v;
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v=t->GetArgs();
      double ke=v[0],th=v[1],vz=v[2],c2=v[5];
      if(c2>chi2Cut||ke<=0) continue;
      if(vz<vzLo||vz>vzHi) continue;
      if(kcOn) th -= (360.0/kcDenom)*(ke-kcPivot);
      KE.push_back(ke); TH_.push_back(th); }
   printf("\n\033[1;33m=== %s : %zu tracks after chi2<%g, vz %g-%g%s ===\033[0m\n",
          label.Data(),KE.size(),chi2Cut,vzLo,vzHi, kcOn?", angle corr ON":"");

   // --- scan Ebeam, anchoring the 6.094 peak -------------------------------------------------
   printf("  Ebeam    mu(6.094)   sigma    Ex(8.317 region)  sigma\n");
   double bestEb=-1,bestD=1e9;
   for(double Eb=156.0;Eb<=168.01;Eb+=0.25){
      TH1D h("h","",380,-2,17);
      for(size_t i=0;i<KE.size();++i){ double e=ExOf(Eb,KE[i],TH_[i]); if(e>-90) h.Fill(e); }
      auto [mu,dmu,sg,ar,ok]=FitPeak(&h,5.55,6.45,0.14);
      if(!ok) continue;
      if(std::fabs(mu-6.094)<bestD){bestD=std::fabs(mu-6.094);bestEb=Eb;}
      if(std::fmod(Eb,1.0)<1e-6){
         auto [m8,d8,s8,a8,ok8]=FitPeak(&h,7.90,9.30,0.25);
         printf("  %6.2f   %7.3f    %5.3f    %s\n",Eb,mu,sg, ok8?Form("%7.3f          %5.3f",m8,s8):"   fit failed");
      }
   }
   if(bestEb<0){printf("\033[1;31mno valid anchor found\033[0m\n");return;}

   // --- at the anchor, measure the 8.317 state ------------------------------------------------
   TH1D *h=new TH1D("hA",Form("^{14}C(p,p') %s, anchored on 6.094;E_{x} [MeV];counts",label.Data()),380,-2,17);
   for(size_t i=0;i<KE.size();++i){ double e=ExOf(bestEb,KE[i],TH_[i]); if(e>-90) h->Fill(e); }
   auto [mu6,d6,s6,a6,ok6]=FitPeak(h,5.55,6.45,0.14);
   auto [mu8,d8,s8,a8,ok8]=FitPeak(h,7.90,9.30,0.25);
   auto [mu0,d0,s0,a0,ok0]=FitPeak(h,-1.2,1.2,0.20);
   printf("\n\033[1;32m  ANCHOR: Ebeam = %.2f MeV\033[0m\n",bestEb);
   printf("   6.094 state -> %.3f +- %.3f  (sigma %.3f, area %.0f) %s\n",mu6,d6,s6,a6,ok6?"":"  <-- FIT SUSPECT");
   printf("   8.317 state -> %.3f +- %.3f  (sigma %.3f, area %.0f) %s\n",mu8,d8,s8,a8,ok8?"":"  <-- FIT SUSPECT");
   printf("   g.s.        -> %.3f +- %.3f  (sigma %.3f, area %.0f) %s\n",mu0,d0,s0,a0,ok0?"":"  <-- FIT SUSPECT");
   if(ok6&&ok8){
      double gain=(mu8-mu6)/(8.317-6.094);
      printf("\n\033[1;33m   residual gain on (Ex - 6.094) = %.4f  (%+.1f %%)\033[0m\n",gain,100*(gain-1));
      printf("   for reference, the old production gave 8.533 i.e. gain 1.098, and the multiplet\n"
             "   spacing gave 1.078 +- 0.015 independently.\n");
   }
   TCanvas *c=new TCanvas("cA","",1000,700); h->GetXaxis()->SetRangeUser(-2,12); h->Draw("hist");
   for(double L:{0.0,6.094,8.317}){auto*l=new TLine(L,0,L,h->GetMaximum()*0.9);l->SetLineColor(kRed+1);l->SetLineStyle(2);l->Draw();}
   c->SaveAs("/home/yassid/a1954_C14_anchor6094.png");
}
