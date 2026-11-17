#include "resfitD.C"
// BOTH pairs tied, each to one amplitude with a fixed internal weight:
//    4.43 + 4.93   (0.50 MeV apart)        5.98 + 6.68   (0.70 MeV apart)
// against a forward-bin response FWHM of 0.93-1.11.  The individual split inside a pair is the
// degenerate direction; the PAIR's yield is what the data determine.
// Check both things at once: fit quality (does the richest bin recover?) and stability (does the
// pair yield survive a 5 % change in the sample?).
void pairs2(){
  setCM(10,5,60);
  RK::exLo=3.0; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true;
  bool keep[6]={true,true,true,true,false,false};
  for(int i=0;i<RK::NR;i++) RK::USE[i]=keep[i];
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    std::vector<std::vector<double>> A[2]; double q[2][16]; double L=0;
    for(int pass=0;pass<2;pass++){
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO = pass? 350 : 300;
      BB_QUIET=true; binbuild();
      A[pass].clear(); double v=rkRun(CFG[ic],NCM,false,&A[pass]);
      if(!pass){ L=v;
        for(int b=0;b<NCM;b++){ RK::Bin& d=RK::B[b]; double s=0; int nd=0;
          for(int i=1;i<=d.h->GetNbinsX();i++){
            double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
            if(x<RK::exLo||x>RK::exHi||n<5) continue;
            double Xt=x-d.gsMu;
            double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                    +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                    +RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
            for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
            if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
          q[0][b]=s/std::max(1,nd-4); } } }
    printf("\n======== %s : both pairs tied ========\n",CFG[ic].nm);
    printf("   -2lnL = %.1f   (%d amplitude slots per bin)\n",L,RK::NGRP);
    printf("   %-10s %10s | %12s %12s | %12s %12s\n",
           "theta_CM","chi2/ndf","D12 yield","refit","D34 yield","refit");
    std::vector<double> r1,r2;
    for(int b=0;b<NCM;b++){
      double a0=(b<(int)A[0].size()&&A[0][b].size()>0)?A[0][b][0]:0;
      double a1=(b<(int)A[1].size()&&A[1][b].size()>0)?A[1][b][0]:0;
      double c0=(b<(int)A[0].size()&&A[0][b].size()>1)?A[0][b][1]:0;
      double c1=(b<(int)A[1].size()&&A[1][b].size()>1)?A[1][b][1]:0;
      printf("   %3.0f-%-6.0f %10.2f | %12.1f %12.1f | %12.1f %12.1f\n",
             TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,q[0][b],a0,a1,c0,c1);
      if(a0>1e-3&&a1>0) r1.push_back(a1/a0);
      if(c0>1e-3&&c1>0) r2.push_back(c1/c0); }
    std::sort(r1.begin(),r1.end()); std::sort(r2.begin(),r2.end());
    if(!r1.empty()) printf("   D12 (4.43+4.93) refit ratio %.2f .. %.2f  (n=%zu/%d)\n",
                           r1.front(),r1.back(),r1.size(),NCM);
    if(!r2.empty()) printf("   D34 (5.98+6.68) refit ratio %.2f .. %.2f  (n=%zu/%d)\n",
                           r2.front(),r2.back(),r2.size(),NCM);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=true; BB_QUIET=true; binbuild();
}
