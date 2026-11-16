#include "resfitC.C"
// The s-wave halo momentum distribution predicts  dsigma/dEx ~ p3* q / Ex^n  with n = EXACTLY 2.
// Float n.  One parameter, same count as the empirical exponential, and the value is physics.
void bkpow(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  for(auto&k:CFG){
    RH::BKG=RH::HALON; double vN=rhRun(k,NCM,false); double n=RH::LASTP0;
    RH::BKG=RH::HALO;  double vH=rhRun(k,NCM,false);
    RH::BKG=RH::EXP;   double vE=rhRun(k,NCM,false);
    printf("\n================ %s : the halo form-factor POWER ================\n",k.nm);
    printf("   p3* q / Ex^n   with n FREE :  n = %.3f        -2lnL = %.1f   (1 par)\n",n,vN);
    printf("   s-wave halo prediction      :  n = 2.000       -2lnL = %.1f   (0 par)\n",vH);
    printf("   empirical exp(-k Ex)        :                  -2lnL = %.1f   (1 par)\n",vE);
    printf("   n - 2 = %+.3f ;  cost of PINNING n at the halo value = %+.1f in -2lnL for 1 dof\n",
           n-2.0, vH-vN);
    printf("   free-n model vs the empirical exponential, SAME parameter count: %+.1f\n",vN-vE);
  }
}
