#include "resfitC.C"
// binbuild.C documents an excess at vz 90-180 mm: ~21 % extra events in that slice, spread
// BROADLY over Ex 1.5-9, i.e. nearly flat.  A flat admixture drags the fitted halo power down.
// Spyral's vz ceiling is 400, so that spike is a far larger fraction of its sample -- which is
// exactly why ATTPCROOT reproduces Spyral's n when restricted to the same window.
// Cut the spike away with VZLO and see whether the two codes converge, and on what.
void bkvzlo(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  double cuts[4]={0,100,200,300};
  printf("\n================ halo power n vs the low-vz cut ================\n");
  printf("   %-12s %10s %10s %10s %10s\n","VZLO [mm]","ATTPCROOT n","-2lnL","Spyral n","-2lnL");
  for(double v:cuts){
    VZLO=v; BB_QUIET=true; binbuild();
    RH::BKG=RH::HALON; double L0=rhRun(CFG[0],NCM,false), n0=RH::LASTP0;
    RH::BKG=RH::HALON; double L1=rhRun(CFG[1],NCM,false), n1=RH::LASTP0;
    printf("   %-12.0f %10.3f %10.1f %10.3f %10.1f\n",v,n0,L0,n1,L1); }
  VZLO=0; BB_QUIET=true; binbuild();
  printf("\n   s-wave halo breakup predicts n = 2 exactly.\n");
}
