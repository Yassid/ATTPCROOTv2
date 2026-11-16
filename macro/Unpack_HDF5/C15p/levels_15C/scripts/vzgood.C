#include "resfitC.C"
// The simulation says the acceptance is FLAT (0.98-1.00) for vz 0-900 and then falls off an
// Ex-DEPENDENT cliff in 900-1000 (0.86 at the g.s. down to 0.35 at Ex = 12).  The Ex calibration
// separately breaks below vz ~ 300.  So the good window is vz in [300,900]: above the calibration
// break, below the acceptance cliff.  ATTPCROOT's adopted 0-1100 includes the cliff; Spyral's
// 0-400 sits in the mis-calibrated region.  Does n stabilise in the good window?
void vzgood(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg b0=CFG[0], b1=CFG[1];
  struct W{const char* nm; double lo,hi;};
  W ws[4]={{"adopted ATTPCROOT 0-1100 / Spyral 0-400",-1,-1},
           {"0-900    (cliff removed only)        ",   0,900},
           {"300-1000 (calib fixed only)          ", 300,1000},
           {"300-900  (BOTH)                      ", 300,900}};
  printf("\n================ continuum exponent n in the good window ================\n");
  printf("   %-40s %16s %16s\n","vz window","ATTPCROOT n","Spyral n");
  for(auto&w:ws){
    if(w.lo<0){ CFG[0]=b0; CFG[1]=b1; VZLO=0; }
    else { CFG[0]=b0; CFG[1]=b1; CFG[0].vzHi=w.hi; CFG[1].vzHi=w.hi; VZLO=w.lo; }
    BB_QUIET=true; binbuild();
    double n[2],s[2];
    for(int ic=0;ic<2;ic++){
      RH::BKG=RH::HALON; double vf=rhRun(CFG[ic],NCM,false); n[ic]=RH::LASTP0;
      RH::SBFIX=true; const double dd=0.20; RH::SBVAL=n[ic]+dd;
      double D=rhRun(CFG[ic],NCM,false)-vf; RH::SBFIX=false;
      s[ic]=(D>1e-6)?dd/std::sqrt(D):0.0; }
    printf("   %-40s %8.3f+-%.3f %8.3f+-%.3f\n",w.nm,n[0],s[0],n[1],s[1]); }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; BB_QUIET=true; binbuild();
  printf("\n   s-wave halo breakup predicts n = 2 exactly.\n");
}
