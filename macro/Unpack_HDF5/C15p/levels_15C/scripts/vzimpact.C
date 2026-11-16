#include "resfitC.C"
// Does removing the vz 90-180 mm excess change the RESONANCE region, or only the continuum?
// Same model both times (parameter-free halo continuum), only the low-vz cut changes.
// If the angular distributions move, the excess was contaminating the states and the cut matters
// for the physics, not just for the background.
void vzimpact(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true; RH::BKG=RH::HALO;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  const double EP[6]={4.43,4.93,5.98,6.68,7.38,8.28};
  for(int ic=0;ic<2;ic++){
    std::vector<std::vector<double>> A0,A2; std::vector<double> c0,c2;
    for(int pass=0;pass<2;pass++){
      VZLO = pass? 200 : 0; BB_QUIET=true; binbuild();
      std::vector<std::vector<double>> out;
      rhRun(CFG[ic],NCM,false,&out);
      std::vector<double> cc;
      for(int b=0;b<NCM;b++){ RH::Bin& d=RH::B[b]; double s=0; int nd=0;
        const double* dummy=nullptr; (void)dummy;
        for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
          if(x<RH::exLo||x>RH::exHi||n<5) continue;
          double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                  +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                  +RH::LASTBG[b]*RH::bshape(x,0.0)*(d.ga?d.ga->Eval(x):1.0);
          for(int r=0;r<RH::NR;r++) if(RH::USE[r]&&d.in[r]) m+=RH::LASTA[b][r]*d.T[r].Eval(x);
          if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
        cc.push_back(s/std::max(1,nd-RH::nper())); }
      if(pass){ A2=out; c2=cc; } else { A0=out; c0=cc; } }
    printf("\n================ %s : effect of VZLO = 200 on the resonances ================\n",CFG[ic].nm);
    printf("   per-angle chi2/ndf   VZLO=0 -> 200 :");
    for(int b=0;b<NCM;b++) printf(" %.2f->%.2f",c0[b],c2[b]); printf("\n");
    double s0=0,s2=0; for(int b=0;b<NCM;b++){s0+=c0[b];s2+=c2[b];}
    printf("   mean chi2/ndf        %.2f -> %.2f\n",s0/NCM,s2/NCM);
    printf("\n   yield ratio  VZLO=200 / VZLO=0, per resonance and theta_CM bin\n");
    printf("   %-8s","Ex");
    for(int b=0;b<NCM;b++) printf(" %5.0f",TCM0+b*TCMSTEP+TCMSTEP/2); printf("\n");
    for(int r=0;r<6;r++){ printf("   %-8.2f",EP[r]);
      for(int b=0;b<NCM;b++){
        double a0=(b<(int)A0.size()&&r<(int)A0[b].size())?A0[b][r]:0;
        double a2=(b<(int)A2.size()&&r<(int)A2[b].size())?A2[b][r]:0;
        if(a0>1e-6) printf(" %5.2f",a2/a0); else printf("     -"); }
      printf("\n"); }
    printf("   (0.83 would be the pure statistics loss; departures from it are real changes)\n");
  }
  VZLO=0; BB_QUIET=true; binbuild();
}
