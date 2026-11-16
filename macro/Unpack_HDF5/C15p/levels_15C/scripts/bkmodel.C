#include "resfitC.C"
// CAN THE CONTINUUM BE MODELLED?
// Every form here is the SOLE baseline (ONETERM), so they genuinely compete instead of sharing
// strength with a separate PS::dens term -- that confound is what made the earlier comparison
// meaningless.  Parameter counts are stated; HALO has ZERO free shape parameters.
void bkmodel(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  struct F{int id; const char* nm; int np;};
  F fs[5]={{RH::PS,   "PS     p3*(M) q(M)                 [uniform 3-body]",0},
           {RH::HALO, "HALO   p3* q / Ex^2                [s-wave halo]   ",0},
           {RH::HALOF,"HALOF  p3* q /(Ex-Sn+Sb)^2         [Sb FLOATED]    ",1},
           {RH::SWAVE,"SWAVE  sqrt(Ex-Sb)/Ex^2            [no phase space]",1},
           {RH::EXP,  "EXP    exp(-k Ex)                  [empirical]     ",1}};
  for(auto&k:CFG){
    printf("\n================ %s : continuum MODEL, single baseline term ================\n",k.nm);
    printf("   %-52s %10s %5s %10s\n","form","-2lnL","npar","shape par");
    double v[5],sp[5]; double best=1e18;
    for(int m=0;m<5;m++){ RH::BKG=fs[m].id; RH::SBFIX=false;
      v[m]=rhRun(k,NCM,false); sp[m]=RH::LASTP0;
      if(v[m]<best) best=v[m]; }
    for(int m=0;m<5;m++)
      printf("   %-52s %10.1f %5d %10.4f  %s\n",fs[m].nm,v[m],fs[m].np,
             fs[m].np?sp[m]:0.0, (v[m]<=best+1e-6)?"<- best":"");
    printf("\n   HALO is parameter-free; Sn(15C) = %.3f MeV.\n",PS::SN);
    printf("   HALOF floats that binding energy -> %.3f MeV  (%+.3f from Sn)\n",
           sp[2],sp[2]-PS::SN);
    printf("   Delta(-2lnL): HALO - best = %+.1f   HALO - EXP = %+.1f (and EXP costs one parameter)\n",
           v[1]-best, v[1]-v[4]);
  }
}
