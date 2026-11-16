#include "resfitC.C"
// 5.98 and 6.68 are 0.70 MeV apart against a forward response FWHM of 0.93-1.11: NOT resolvable,
// and the anatomy shows them trading (at theta_CM 15-20 the 5.98 tail gives 14.32 at Ex = 6.68
// while the 6.68 itself fits to 0.00).  That is the forward collapse.
// Does dropping 5.98 recover a stable 6.68 across ALL ten angular bins?
void pair598(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALON; RH::TIE12=true;
  Cfg b0=CFG[0], b1=CFG[1];
  struct S{const char* nm; bool u[6];};
  S sets[2]={{"D12 + 5.98 + 6.68 (as before)",{true,true,true,true,false,false}},
             {"D12 + 6.68 only              ",{true,true,false,true,false,false}}};
  for(auto&st:sets){
    printf("\n============ resonance set: %s ============\n",st.nm);
    for(int ic=0;ic<2;ic++){
      std::vector<std::vector<double>> A[2];
      for(int pass=0;pass<2;pass++){
        for(int i=0;i<RH::NR;i++) RH::USE[i]=st.u[i];
        CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO = pass? 350 : 300;
        BB_QUIET=true; binbuild();
        A[pass].clear(); rhRun(CFG[ic],NCM,false,&A[pass]); }
      // 6.68 is the LAST entry of the row in both sets
      printf("   %-10s  6.68 yield per theta_CM bin (and refit at VZLO 350)\n",CFG[ic].nm);
      printf("      %-10s","theta_CM"); for(int b=0;b<NCM;b++) printf(" %6.0f",TCM0+(b+0.5)*TCMSTEP);
      printf("\n      %-10s","yield");
      int nz=0,ng=0;
      std::vector<double> rat;
      for(int b=0;b<NCM;b++){
        int k=(b<(int)A[0].size())?(int)A[0][b].size()-1:-1;
        double y0=(k>=0)?A[0][b][k]:0;
        int k1=(b<(int)A[1].size())?(int)A[1][b].size()-1:-1;
        double y1=(k1>=0)?A[1][b][k1]:0;
        printf(" %6.1f",y0); if(y0<1e-3) nz++; else { ng++; if(y1>0) rat.push_back(y1/y0); } }
      printf("\n      %-10s","refit");
      for(int b=0;b<NCM;b++){ int k1=(b<(int)A[1].size())?(int)A[1][b].size()-1:-1;
        printf(" %6.1f",(k1>=0)?A[1][b][k1]:0); }
      std::sort(rat.begin(),rat.end());
      printf("\n      bins with yield: %d/%d   zeros: %d", ng, NCM, nz);
      if(!rat.empty()) printf("   refit ratio %.2f .. %.2f (median %.2f)",
                              rat.front(),rat.back(),rat[rat.size()/2]);
      printf("\n");
    }
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true; BB_QUIET=true; binbuild();
}
