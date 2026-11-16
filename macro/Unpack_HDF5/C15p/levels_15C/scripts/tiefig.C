#include "resfitC.C"
// STABILITY of the angular distributions, which ARE the physics deliverable.
// Refitting the same data under a small change (here the vz window) moved individual per-bin
// yields by factors from 0.00 to 94467 -- an ill-conditioned fit, not physics.  The 4.43/4.93
// pair sits 0.5 MeV apart against a response FWHM of ~0.93, so it is NOT resolvable and the two
// amplitudes trade freely.  TIE12 replaces them by ONE doublet amplitude with a shared internal
// ratio.  Does that make the yields reproducible?
//   metric: for each resonance, the spread of (yield at window B)/(yield at window A) across the
//   theta_CM bins.  A stable fit gives a tight cluster near the statistics ratio.
void tiefig(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALO;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg b0=CFG[0], b1=CFG[1];
  const double EP[6]={4.43,4.93,5.98,6.68,7.38,8.28};
  for(int tie=0;tie<2;tie++){
    RH::TIE12 = (tie==1);
    printf("\n================ angular-distribution stability : TIE12 %s ================\n",
           tie? "ON":"OFF");
    for(int ic=0;ic<2;ic++){
      std::vector<std::vector<double>> A[2];
      for(int pass=0;pass<2;pass++){
        CFG[0]=b0; CFG[1]=b1;
        CFG[ic].vzHi = 900; VZLO = pass? 350 : 300;     // a DELIBERATELY SMALL change
        BB_QUIET=true; binbuild();
        A[pass].clear(); rhRun(CFG[ic],NCM,false,&A[pass]); }
      printf("   %-10s  ratio of yields, VZLO 350 vs 300 (a ~5%% change in statistics)\n",CFG[ic].nm);
      int nr = RH::TIE12 ? 5 : 6;
      for(int r=0;r<nr;r++){
        std::vector<double> v;
        for(int b=0;b<NCM;b++){
          double a0=(b<(int)A[0].size()&&r<(int)A[0][b].size())?A[0][b][r]:0;
          double a1=(b<(int)A[1].size()&&r<(int)A[1][b].size())?A[1][b][r]:0;
          if(a0>1e-3) v.push_back(a1/a0); }
        if(v.empty()){ printf("      %-6s   (no yield)\n", RH::TIE12&&r==0?"D12":Form("%.2f",EP[RH::TIE12?r+1:r])); continue; }
        std::sort(v.begin(),v.end());
        double med=v[v.size()/2], lo=v.front(), hi=v.back();
        printf("      %-6s  n=%2zu   median %6.2f   range %7.3f .. %9.3f\n",
               RH::TIE12&&r==0?"D12":Form("%.2f",EP[RH::TIE12?r+1:r]),v.size(),med,lo,hi); }
    }
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false; BB_QUIET=true; binbuild();
  printf("\n   A stable fit gives a tight cluster near 1; a wide range means the amplitude is\n"
         "   not determined by the data and the angular distribution is not a measurement.\n");
}
