#include "resfitD.C"
// Pin the seventh state, and test whether it is an artefact of the fit window's lower edge:
// a component sitting near Ex = exLo can absorb edge effects.  Move exLo down to 2.5 and 2.0;
// a real state keeps its energy, an edge artefact follows the edge.
void ex37(){
  setCM(10,5,60);
  RK::ONETERM=true; RK::BKG=RK::HALON; RK::TIE12=false; RK::PAIRS=true;
  Cfg b0=CFG[0], b1=CFG[1];
  double lows[3]={3.0,2.5,2.0};
  for(int ic=0;ic<2;ic++){
    printf("\n======== %s : where is the seventh state? ========\n",CFG[ic].nm);
    printf("   %-12s %10s %10s %12s\n","fit from Ex","best Ex","gain","(step 0.05)");
    for(double lo:lows){
      RK::exLo=lo; RK::exHi=9.0;
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
      bool off[7]={true,true,true,true,false,false,false};
      for(int i=0;i<RK::NR;i++) RK::USE[i]=off[i];
      double L0=rkRun(CFG[ic],NCM,false);
      double best=1e18,bE=0;
      for(double E=3.30;E<=4.21;E+=0.05){
        RK::EPP6[6]=E+RK::PPGS6;
        bool on[7]={true,true,true,true,false,false,true};
        for(int i=0;i<RK::NR;i++) RK::USE[i]=on[i];
        double L=rkRun(CFG[ic],NCM,false);
        if(L<best){best=L;bE=E;} }
      printf("   %-12.1f %10.2f %10.1f\n",lo,bE,L0-best); }
  }
  RK::exLo=3.0; CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6); BB_QUIET=true; binbuild();
}
