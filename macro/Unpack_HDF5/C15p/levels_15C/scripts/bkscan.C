#include "resfitC.C"
// 1) Profile -2lnL in the halo power n to get a real uncertainty (Delta = 1 is 1 sigma).
// 2) Does the PHYSICAL continuum improve the per-angle fits, which is the standing complaint?
void bkscan(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  for(auto&k:CFG){
    printf("\n================ %s : profile of the halo power n ================\n",k.nm);
    RH::BKG=RH::HALON; double vfree=rhRun(k,NCM,false), nhat=RH::LASTP0;
    printf("   n_hat = %.3f   -2lnL = %.1f\n",nhat,vfree);
    printf("   %6s %10s %8s\n","n","-2lnL","Delta");
    double lo=-9,hi=-9;
    for(double n=0.6;n<=3.01;n+=0.2){
      RH::BKG=RH::HALON; RH::SBFIX=true; RH::SBVAL=n;   // SBFIX path fixes p[0]
      double v=rhRun(k,NCM,false); RH::SBFIX=false;
      printf("   %6.2f %10.1f %8.1f %s\n",n,v,v-vfree,(std::fabs(v-vfree-1.0)<0.35)?"<- 1 sigma":"");
      if(v-vfree<1.0){ if(lo<-8) lo=n; hi=n; } }
    printf("   => n = %.2f  (1-sigma band roughly %.2f .. %.2f)   s-wave halo predicts 2\n",nhat,lo,hi);
  }
  // per-angle quality, physical model vs empirical
  for(auto&k:CFG){
    printf("\n================ %s : per-angle chi2/ndf ================\n",k.nm);
    RH::BKG=RH::EXP;   printf("   --- empirical exp(-k Ex) ---\n");   rhRun(k,NCM,true);
    RH::BKG=RH::HALON; printf("\n   --- physical p3* q / Ex^n ---\n"); rhRun(k,NCM,true);
  }
}
