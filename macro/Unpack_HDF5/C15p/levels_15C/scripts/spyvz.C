#include "resfitC.C"
// Spyral's adopted vzHi = 400 cuts at its own median (508): more than half the events, the upper
// half, are discarded.  It was adopted as a workaround for the Ebeam(vz) knob.  Is it still
// needed?  Widen it and watch the two things that would break first:
//   - the g.s. CENTROID (the calibration diagnostic -- it must stay at 0)
//   - the g.s. WIDTH (a wrong Ebeam(vz) slope smears the peak as the window widens)
// and then see what the continuum power does with NO low-vz cut at all.
void spyvz(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg base=CFG[1];
  double hi[4]={400,600,800,1000};
  printf("\n================ Spyral : widening the vz window (VZLO = 0) ================\n");
  printf("   %-8s %9s   %9s %9s   %9s %9s\n","vzHi","N kept","gs centroid","gs FWHM","halo n","-2lnL");
  for(double h:hi){
    CFG[1]=base; CFG[1].vzHi=h; VZLO=0; BB_QUIET=true; binbuild();
    // g.s. peak on the full spectrum of the richest bin
    TFile f("dd_cmbins.root"); double N=0; TH1D* hb=nullptr; double bN=-1;
    for(int b=0;b<NCM;b++){ TH1D* x=(TH1D*)f.Get(Form("Spyral_cm%d",b));
      if(!x) continue; N+=x->GetEntries(); if(x->GetEntries()>bN){bN=x->GetEntries(); hb=(TH1D*)x->Clone("hb"); hb->SetDirectory(nullptr);} }
    f.Close();
    double mu=0,fw=0;
    if(hb){ TF1 g("g","gaus",-1.2,1.2); hb->Fit(&g,"RQ0N");
            mu=g.GetParameter(1); fw=2.355*g.GetParameter(2); }
    RH::BKG=RH::HALON; double L=rhRun(CFG[1],NCM,false);
    printf("   %-8.0f %9.0f   %9.3f %9.3f   %9.3f %9.1f\n",h,N,mu,fw,RH::LASTP0,L); }
  CFG[1]=base; VZLO=0; BB_QUIET=true; binbuild();
  printf("\n   ATTPCROOT, for reference, uses vzHi = 1100 and returns n = 1.832 at VZLO = 0.\n");
}
