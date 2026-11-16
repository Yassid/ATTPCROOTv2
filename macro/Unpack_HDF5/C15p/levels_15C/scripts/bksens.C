#include "resfitC.C"
// Why do the two codes return different halo powers (1.83 vs 1.22)?  Their SELECTIONS differ a
// lot -- chi2 cut 0.2 vs 5.0, vz ceiling 1100 vs 400, Ebeam 191.25 vs 195 -- and all three enter
// n, through the admitted track quality, through the acceptance and through PS::dens.
// Move ATTPCROOT's cuts toward Spyral's one at a time and watch n.
void bksens(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg base=CFG[0];
  struct V{const char* nm; double chi2,vzHi,eb;};
  V vs[6]={{"ATTPCROOT as adopted        ",0.2 ,1100,191.25},
           {"  chi2 0.2 -> 0.5           ",0.5 ,1100,191.25},
           {"  chi2 0.2 -> 1.0           ",1.0 ,1100,191.25},
           {"  vz ceiling 1100 -> 400    ",0.2 , 400,191.25},
           {"  Ebeam 191.25 -> 195       ",0.2 ,1100,195.00},
           {"  all three, Spyral-like    ",1.0 , 400,195.00}};
  printf("\n================ sensitivity of the halo power n to the selection ================\n");
  printf("   %-30s %7s %7s %8s   %8s %9s\n","variant","chi2","vzHi","Ebeam","n","-2lnL");
  for(auto&v:vs){
    CFG[0]=base; CFG[0].chi2=v.chi2; CFG[0].vzHi=v.vzHi; CFG[0].eb=v.eb;
    BB_QUIET=true; binbuild();            // the binned spectra MUST be regenerated for new cuts
    RH::BKG=RH::HALON; double L=rhRun(CFG[0],NCM,false);
    printf("   %-30s %7.1f %7.0f %8.2f   %8.3f %9.1f\n",v.nm,v.chi2,v.vzHi,v.eb,RH::LASTP0,L); }
  CFG[0]=base; BB_QUIET=true; binbuild();  // restore
  printf("\n   Spyral returns n = 1.220 +- 0.07 on its own selection; ATTPCROOT 1.832 +- 0.10.\n");
}
