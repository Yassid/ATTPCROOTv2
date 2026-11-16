#include "resfitC.C"
// 7.38 and 8.28 are not determined by the data: their yields swing through 0.000 under a 5 %
// change in statistics, in BOTH codes, with or without TIE12.  Carrying an undetermined component
// lets it absorb whatever the continuum misses, which can destabilise the ones that ARE measured.
// Does dropping them stabilise D12, 5.98 and 6.68?
void redset(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALO; RH::TIE12=true;
  Cfg b0=CFG[0], b1=CFG[1];
  const char* nm[5]={"D12","5.98","6.68","7.38","8.28"};
  for(int set=0;set<2;set++){
    bool keep[6]={true,true,true,true, set==0, set==0};   // set 1 drops 7.38 and 8.28
    printf("\n============ TIE12 ON, resonance set = %s ============\n",
           set? "D12 + 5.98 + 6.68 only":"all six");
    for(int ic=0;ic<2;ic++){
      std::vector<std::vector<double>> A[2];
      for(int pass=0;pass<2;pass++){
        for(int i=0;i<RH::NR;i++) RH::USE[i]=keep[i];
        CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO = pass? 350 : 300;
        BB_QUIET=true; binbuild();
        A[pass].clear(); rhRun(CFG[ic],NCM,false,&A[pass]); }
      printf("   %-10s  yield ratio VZLO 350/300, per resonance\n",CFG[ic].nm);
      int nr = set? 3 : 5;
      for(int r=0;r<nr;r++){
        std::vector<double> v;
        for(int b=0;b<NCM;b++){
          double a0=(b<(int)A[0].size()&&r<(int)A[0][b].size())?A[0][b][r]:0;
          double a1=(b<(int)A[1].size()&&r<(int)A[1][b].size())?A[1][b][r]:0;
          if(a0>1e-3) v.push_back(a1/a0); }
        if(v.empty()){ printf("      %-6s (no yield)\n",nm[r]); continue; }
        std::sort(v.begin(),v.end());
        int nz=0; for(double x:v) if(x<1e-6) nz++;
        printf("      %-6s n=%2zu  median %5.2f  range %6.3f .. %8.3f   zeros: %d\n",
               nm[r],v.size(),v[v.size()/2],v.front(),v.back(),nz); }
    }
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true; BB_QUIET=true; binbuild();
}
