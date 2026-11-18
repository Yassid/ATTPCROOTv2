#include "resfitD.C"
// The disputed "state at ~3.5" sits at the LOW EDGE of the fit window, exactly where a component
// absorbs any model deficiency.  Rather than argue about it: cut the window above it.
// For each lower limit, report fit quality, whether a 7th component is still wanted, and whether
// the D34 (5.98+6.68) yield -- the deliverable -- survives a 5 % change in the sample.
void lowcut(){
  setCM(10,5,60);
  RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON; RK::TIE12=false; RK::PAIRS=true;
  RK::NPERBIN=false;
  Cfg b0=CFG[0], b1=CFG[1];
  double los[4]={3.0,3.5,4.0,4.2};
  for(int ic=0;ic<2;ic++){
    printf("\n================ %s : where should the fit start? ================\n",CFG[ic].nm);
    printf("   %-9s %9s %9s %9s   %-22s %s\n",
           "Ex from","mean c2","c2 15-20","worst c2","7th state gain","D34 refit spread");
    for(double lo:los){
      RK::exLo=lo;
      bool base[7]={true,true,true,true,false,false,false};
      std::vector<std::vector<double>> A[2];
      double mc=0,c1520=0,worst=0; int nb=0;
      for(int pass=0;pass<2;pass++){
        for(int i=0;i<RK::NR;i++) RK::USE[i]=base[i];
        CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO = pass? 350 : 300;
        BB_QUIET=true; binbuild();
        A[pass].clear(); double L=rkRun(CFG[ic],NCM,false,&A[pass]);
        if(pass) continue;
        for(int b=0;b<NCM;b++){ RK::Bin& d=RK::B[b]; if(!d.h) continue;
          double s=0; int nd=0;
          for(int i=1;i<=d.h->GetNbinsX();i++){
            double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
            if(x<RK::exLo||x>RK::exHi||n<5) continue;
            double Xt=x-d.gsMu;
            double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                    +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                    +RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
            for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
            if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
          if(nd<6) continue;
          double q=s/std::max(1,nd-4); mc+=q; nb++;
          if(b==1) c1520=q; worst=std::max(worst,q); } }
      // is a 7th component still wanted?
      for(int i=0;i<RK::NR;i++) RK::USE[i]=base[i];
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
      double L0=rkRun(CFG[ic],NCM,false);
      double best=1e18;
      for(double E=3.20;E<=4.21;E+=0.15){
        if(E<RK::exLo-0.4) continue;
        RK::EPP6[6]=E+RK::PPGS6;
        bool on[7]={true,true,true,true,false,false,true};
        for(int i=0;i<RK::NR;i++) RK::USE[i]=on[i];
        double L=rkRun(CFG[ic],NCM,false); if(L<best) best=L; }
      double gain = (best<1e17)? L0-best : 0;
      // D34 stability
      std::vector<double> r;
      for(int b=0;b<NCM;b++){
        double a0=(b<(int)A[0].size()&&A[0][b].size()>1)?A[0][b][1]:0;
        double a1=(b<(int)A[1].size()&&A[1][b].size()>1)?A[1][b][1]:0;
        if(a0>1e-3&&a1>0) r.push_back(a1/a0); }
      std::sort(r.begin(),r.end());
      printf("   %-9.1f %9.2f %9.2f %9.2f   %-22s %s\n",lo,nb?mc/nb:0,c1520,worst,
             Form("%+.1f  (+10 par)",gain),
             r.empty()?"-":Form("%.2f .. %.2f  (n=%zu/%d)",r.front(),r.back(),r.size(),NCM)); }
  }
  RK::exLo=3.0; CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6); BB_QUIET=true; binbuild();
}
