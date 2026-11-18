#include "resfitD.C"
// A DIFFERENT BACKGROUND SHAPE PER ANGULAR BIN, states fixed.
// The continuum is the 1n breakup channel and its Ex dependence is the inelastic form factor,
// which depends on momentum transfer -- and q varies strongly across theta_CM 10-60.  Forcing one
// shape on every angle is wrong, and shows up as chi2 elevated UNIFORMLY rather than failing in
// one place.
// The first attempt at this returned a WORSE likelihood for a strictly larger model: a cold start
// with 10 extra free exponents never converged.  Fix: fit the global-shape model first, then
// release the per-bin exponents SEEDED from that solution.
void perbin(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  auto chis=[&](double* out){
    double mc=0; int nb=0;
    for(int b=0;b<NCM;b++){ RK::Bin& d=RK::B[b]; if(!d.h){ out[b]=0; continue; }
      double s=0; int nd=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){
        double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
        if(x<RK::exLo||x>RK::exHi) continue;
        double Xt=x-d.gsMu;
        double shp = RK::NPERBIN ? RK::LASTP[1+b*RK::nper()+2+RK::NGRP] : RK::LASTP0;
        double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                +RK::LASTBG[b]*RK::bshape(Xt,shp)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
        for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
        if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/std::max(m,1.0); nd++; }
      out[b]=s/std::max(1.0,nd-5.3); mc+=out[b]; nb++; }
    return nb? mc/nb : 0; };
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    // --- pass 1: one global shape
    RK::NPERBIN=false; RK::SEED.clear();
    double L1=rkRun(CFG[ic],NCM,false);
    double c1[16]; double m1=chis(c1);
    std::vector<double> sol=RK::LASTP; double nglob=RK::LASTP0;
    int nper0=RK::nper();
    // --- pass 2: per-bin shape, seeded from pass 1
    RK::NPERBIN=true;
    int nper1=RK::nper();                       // one more slot per bin
    std::vector<double> sd(1+NCM*nper1+RK::NRAT,0.0);
    sd[0]=0.0;
    for(int b=0;b<NCM;b++){
      for(int j=0;j<nper0;j++) sd[1+b*nper1+j]=sol[1+b*nper0+j];
      sd[1+b*nper1+nper0]=nglob; }              // each bin starts at the GLOBAL exponent
    for(int r=0;r<RK::NRAT;r++) sd[1+NCM*nper1+r]=sol[1+NCM*nper0+r];
    RK::SEED=sd;
    double L2=rkRun(CFG[ic],NCM,false);
    double c2[16]; double m2=chis(c2);
    RK::SEED.clear(); RK::NPERBIN=false;
    printf("\n================ %s : one background shape, or one per angle? ================\n",
           CFG[ic].nm);
    printf("   global shape   -2lnL %9.1f   mean chi2 %.2f   (n = %.2f)\n",L1,m1,nglob);
    printf("   per-bin shape  -2lnL %9.1f   mean chi2 %.2f   gain %+.1f for %d parameters\n",
           L2,m2,L1-L2,NCM);
    printf("\n   %-10s %10s %10s %10s\n","theta_CM","chi2 glob","chi2 /bin","fitted n");
    for(int b=0;b<NCM;b++)
      printf("   %3.0f-%-6.0f %10.2f %10.2f %10.2f\n",
             TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,c1[b],c2[b],
             sd.size()? RK::LASTP[1+b*nper1+nper0] : 0.0);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
