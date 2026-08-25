/// @file multiplet_gain_C14.C
/// @brief Second, independent determination of the 14C excitation scale from the SPACING of the
///        6-7.5 MeV multiplet: member energies pinned, one common gain free.
///
/// WIDTH MODEL. 12_excited.tex fixes every component at one width, because floated freely they run
/// to 0.20-0.42 MeV and the amplitudes go anti-correlated. But a single fixed width is also wrong:
/// the resolution DEGRADES with excitation energy (the same effect that makes a fixed Ex window
/// distort the g.s. yield, 09_extraction.tex), so the upper members are genuinely broader than the
/// 6.094 anchor. Fixing all four at sigma(6.094) forces the fit to absorb that in the amplitudes
/// and the background -- which is what gave chi2/ndf 3.86 on the 2026-08-25 CATIMA production.
///
/// The compromise used here: sigma(E) = s0 + s1*(E - 6.094), TWO free parameters for four
/// components. The width still grows with energy, but it cannot be tuned per peak, so the
/// anti-correlation that forced the original fixing cannot come back.
///
///   root -b -q 'multiplet_gain_C14.C("plots/proton_kin_cmpB_cat.root",159.75)'
static const double UU = 931.49401;
static double om2b(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
static double ExOfB(double Eb,double K,double thDeg)
{
   const double m1=14.003242*UU,m2=1.007825*UU,m3=1.007825*UU,m4=14.003242*UU;
   double Et1=Eb+m1,Et3=K+m3,s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double th=thDeg*TMath::DegToRad();
   double a=(std::cos(th)*om2b(s,m1*m1,m2*m2)*om2b(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   return a<0?-99.0:std::sqrt(a)-m4;
}
static double gE0[4]={6.09,6.70,7.00,7.27};
static int    gNbg = 1;      // 1 = linear background, 2 = quadratic
static bool   gVarW = true;  // width grows with Ex
// p: [0]=gain [1..4]=amplitudes [5]=s0 [6]=s1 [7..]=background
static double multiF(double *x,double *p)
{
   double s = p[7] + p[8]*x[0] + (gNbg>=2 ? p[9]*x[0]*x[0] : 0.0);
   for (int i=0;i<4;i++){
      double mu  = 6.094 + p[0]*(gE0[i]-6.094);
      double sig = p[5] + (gVarW ? p[6]*(gE0[i]-6.094) : 0.0);
      if (sig < 0.02) sig = 0.02;
      s += p[1+i]*std::exp(-0.5*std::pow((x[0]-mu)/sig,2));
   }
   return s;
}
void multiplet_gain_C14(TString cache="plots/proton_kin_cmpB_cat.root", double Ebeam=159.75,
                        double sig0Seed=0.111, double chi2Cut=5.0, double vzLo=20, double vzHi=1100,
                        double kcDenom=6428.6, double kcPivot=0.0)
{
   gStyle->SetOptStat(0);
   TFile *f=TFile::Open(cache); TNtuple *t=(TNtuple*)f->Get("pk"); float *v;
   TH1D *h=new TH1D("hM","^{14}C multiplet, common-scale fit;E_{x} [MeV];counts",300,-2,16);
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v=t->GetArgs();
      if(v[5]>chi2Cut||v[0]<=0||v[2]<vzLo||v[2]>vzHi) continue;
      double th=v[1]-(360.0/kcDenom)*(v[0]-kcPivot);
      double e=ExOfB(Ebeam,v[0],th); if(e>-90) h->Fill(e); }

   for (int варW=1; варW>=0; --варW) {          // varying width, then the old fixed-width control
      for (int nbg=1; nbg<=2; ++nbg) {
         gVarW = варW; gNbg = nbg;
         TF1 F("F",multiF,5.60,7.70,10);
         F.SetParameters(1.0, 80,250,190,110, sig0Seed, 0.05, 5,0,0);
         F.SetParLimits(0,0.90,1.25);
         for(int i=1;i<=4;i++) F.SetParLimits(i,0,2000);
         F.SetParLimits(5,0.05,0.45);
         if (варW) F.SetParLimits(6,-0.10,0.40); else F.FixParameter(6,0.0);
         if (nbg<2) F.FixParameter(9,0.0);
         TFitResultPtr r=h->Fit(&F,"RQNS");
         double c2n = (r.Get()&&r->Ndf()>0) ? r->Chi2()/r->Ndf() : -1;
         printf("\n\033[1;33m width %s, background %s : chi2/ndf = %.2f\033[0m\n",
                варW?"GROWS with Ex":"FIXED (the old model)", nbg==1?"linear":"quadratic", c2n);
         printf("   gain   = %.4f +- %.4f%s\n", F.GetParameter(0), F.GetParError(0),
                (std::fabs(F.GetParameter(0)-0.90)<1e-3||std::fabs(F.GetParameter(0)-1.25)<1e-3)?"  <-- RAILED":"");
         printf("   sigma0 = %.4f +- %.4f   dsigma/dEx = %+.4f +- %.4f%s\n",
                F.GetParameter(5),F.GetParError(5),F.GetParameter(6),F.GetParError(6),
                (варW&&(std::fabs(F.GetParameter(6)-(-0.10))<1e-4||std::fabs(F.GetParameter(6)-0.40)<1e-4))?"  <-- RAILED":"");
         for(int i=0;i<4;i++)
            printf("     %.2f -> %.3f  sigma %.3f  A = %6.1f +- %5.1f%s\n", gE0[i],
                   6.094+F.GetParameter(0)*(gE0[i]-6.094),
                   F.GetParameter(5)+(варW?F.GetParameter(6)*(gE0[i]-6.094):0.0),
                   F.GetParameter(1+i),F.GetParError(1+i),
                   F.GetParameter(1+i)<1e-3?"  <-- COLLAPSED":"");
         if (варW && nbg==1) {
            TCanvas *c=new TCanvas("cM","",1000,650);
            h->GetXaxis()->SetRangeUser(4.8,8.6); h->Draw("hist");
            F.SetLineColor(kRed+1); F.SetNpx(700); F.DrawCopy("same");
            c->SaveAs("/home/yassid/a1954_C14_multiplet_gain.png");
         }
      }
   }
}
