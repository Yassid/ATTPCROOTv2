#include "resfitC.C"
// THE CONTINUUM MODEL, with the vz 90-180 mm instrumental excess removed (VZLO = 200 mm).
//   direct s-wave halo breakup:   dsigma/dEx  ~  p3*(M) q(M) / Ex^2     ZERO free parameters
// Profile the power for an uncertainty, then put the parameter-free model head to head with the
// empirical exponential, which costs one parameter.
void bkfinal(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  VZLO=200; BB_QUIET=true; binbuild();
  for(auto&k:CFG){
    printf("\n================ %s   VZLO = 200 mm ================\n",k.nm);
    RH::BKG=RH::HALON; double vf=rhRun(k,NCM,false), nh=RH::LASTP0;
    // 1-sigma from the profile
    double s=0; { RH::SBFIX=true;
      for(double d=0.02;d<1.2;d+=0.02){ RH::SBVAL=nh+d; double v=rhRun(k,NCM,false);
        if(v-vf>=1.0){ s=d; break; } } RH::SBFIX=false; }
    RH::BKG=RH::HALO; double vH=rhRun(k,NCM,false);     // n = 2 exactly, 0 par
    RH::BKG=RH::EXP;  double vE=rhRun(k,NCM,false);     // empirical, 1 par
    RH::BKG=RH::PS;   double vP=rhRun(k,NCM,false);     // uniform 3-body, 0 par
    printf("   p3* q / Ex^n,  n FREE        n = %.3f +- %.3f      -2lnL = %8.1f   (1 par)\n",nh,s,vf);
    printf("   s-wave halo PREDICTION       n = 2 exactly          -2lnL = %8.1f   (0 par)\n",vH);
    printf("   empirical exp(-kappa Ex)                            -2lnL = %8.1f   (1 par)\n",vE);
    printf("   uniform 3-body p3* q                                -2lnL = %8.1f   (0 par)\n",vP);
    printf("\n   n - 2 = %+.3f  = %.1f sigma\n",nh-2.0,std::fabs(nh-2.0)/std::max(1e-6,s));
    printf("   parameter-free halo model vs the empirical exponential: %+.1f in -2lnL,\n"
           "   with ONE FEWER parameter.\n",vH-vE);
    printf("   halo form factor buys %.1f over uniform phase space for NOTHING.\n",vP-vH);
  }
  VZLO=0; BB_QUIET=true; binbuild();
}
