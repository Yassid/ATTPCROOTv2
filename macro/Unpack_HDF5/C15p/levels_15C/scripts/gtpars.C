#include "resfitD.C"
// What IS the fitted response, bin by bin?  alpha is where the Gaussian core hands over to the
// exponential tail, in units of sigma.  alpha ~ 1 means there is essentially NO Gaussian core --
// the "response" is a two-sided exponential.  Print it, with the g.s. amplitude, so the claim
// "the tails only matter where the g.s. is large" can be checked rather than asserted.
void gtpars(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::B.clear(); RK::refValid=false; RK::refN=0;
    for(int b=0;b<NCM;b++) rkPrep(CFG[ic],b);
    RK::B.clear();
    for(int b=0;b<NCM;b++) rkPrep(CFG[ic],b);
    printf("\n======== %s : the fitted GT response, per angle ========\n",CFG[ic].nm);
    printf("   %-9s %9s %8s %8s %8s %8s %9s %9s\n",
           "theta_CM","gs amp","FWHM_L","FWHM_R","alpha_L","alpha_R","tail@Ex4","tail@Ex6");
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      double t4=d.gsAmp*shapeGT(4.0,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
      double t6=d.gsAmp*shapeGT(6.0,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
      printf("   %3.0f-%-5.0f %9.0f %8.2f %8.2f %8.2f %8.2f %9.2f %9.2f\n",
             d.cmLo,d.cmHi,d.gsAmp,2.355*d.sgL,2.355*d.sgR,d.alL,d.alR,t4,t6); }
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
