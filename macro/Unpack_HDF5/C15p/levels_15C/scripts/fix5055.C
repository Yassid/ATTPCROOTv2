#include "angfinal.C"
// Spyral 50-55 is the one genuine failure: chi2R 920, EDM 2e+05, status 303, with the background
// collapsing vertically near Ex ~ 2.4 while the states spike to compensate.  What is different
// about that bin?  Print its spectrum and what the fit is doing there.
void fix5055(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0; B1::exLo=-1.5;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  int ic=1, b=8;                       // Spyral, theta_CM 50-55
  CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900;
  VZLO=300; BB_QUIET=true; binbuild();
  RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
  for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
  RK::Bin& d=RK::B[b];
  printf("\n======== Spyral theta_CM %.0f-%.0f ========\n",d.cmLo,d.cmHi);
  printf("   N = %.0f   gsAmp = %.0f   gsMu = %+.3f   FWHM_R = %.2f\n",
         d.N,d.gsAmp,d.gsMu,2.355*d.sgR);
  printf("   templates in window:"); for(int r=0;r<RK::NR;r++)
    printf(" %.2f=%s",RK::EPP6[r]-RK::PPGS6,d.in[r]?"Y":"n"); printf("\n");
  printf("\n   %6s %8s   %6s %8s\n","Ex","counts","Ex","counts");
  for(double X=-1.5;X<=4.0;X+=0.25){
    int i1=d.h->FindBin(X), i2=d.h->FindBin(X+5.75);
    printf("   %6.2f %8.0f   %6.2f %8.0f\n",X,d.h->GetBinContent(i1),
           X+5.75,d.h->GetBinContent(i2)); }
  double tot=0; for(int i=1;i<=d.h->GetNbinsX();i++){
    double x=d.h->GetBinCenter(i); if(x>=-1.5&&x<=10) tot+=d.h->GetBinContent(i); }
  printf("\n   total counts in the fit window: %.0f over %d bins\n",tot,
         (int)((10.0+1.5)/0.1));
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
