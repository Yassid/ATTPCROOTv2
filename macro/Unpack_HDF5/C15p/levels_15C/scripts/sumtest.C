#include "simfit.C"
// Is the TOTAL inelastic strength per theta_CM bin smooth even when the per-resonance
// decomposition is not?  If so, the data measures the SUM and not the decomposition -- which
// would explain the bin-to-bin jumps, the background systematic, and the free-energy migration
// all at once.
//
// Integrate each fitted model over the resonance region (Ex 3.5-9) with the continuum and
// background subtracted, and compare against the same integral taken directly from the DATA.
void sumtest(){
  for(auto&k:CFG){
    double nll=sfRun(k,false,false);
    const int NB2=SF::NBIN;
    printf("\n== %s  (free per-bin amplitudes, -2lnL %.1f)\n",k.nm,nll);
    printf("   bin  cm      N  | sum of BW areas | data minus (cont+bkg+g.s.) | ratio\n");
    // recover the parameters by re-running and caching: sfRun already left SF::B prepared,
    // so refit once more and keep the minimiser result
    SF::SMOOTH=false;
    int nper=SF::NPB, npar=2+NB2*nper;
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000);
    m->SetTolerance(0.1); m->SetPrintLevel(0);
    ROOT::Math::Functor fnc(&SF::fcn,npar); m->SetFunction(fnc);
    m->SetLimitedVariable(0,"bkg_w",1.0,0.02,0.05,3.0);
    m->SetLimitedVariable(1,"bkg_kappa",0.3,0.01,0.0,2.0);
    for(int b=0;b<NB2;b++){ SF::BinData&d=SF::B[b]; if(!d.h) continue;
      double pk=d.h->GetBinContent(d.h->GetMaximumBin());
      double lvl=0;int nl=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){double x=d.h->GetBinCenter(i);
        if(x>6&&x<8){lvl+=d.h->GetBinContent(i);nl++;}} lvl=nl?lvl/nl:1;
      int o=2+b*nper;
      m->SetLimitedVariable(o+0,Form("gs_%d",b),pk,pk*1e-3,0,pk*3);
      m->SetLimitedVariable(o+1,Form("e52_%d",b),pk*0.05,pk*1e-4,0,pk);
      m->SetLimitedVariable(o+2,Form("cont_%d",b),1e-3,1e-5,0,1e3);
      m->SetLimitedVariable(o+3,Form("bkgB_%d",b),0.5*lvl,lvl*1e-3,0,lvl*20);
      for(int i=0;i<6;i++){ double sc=d.in[i]?d.T[i].Eval((EPP[i]-PPGS)+d.gsMu):0;
        double seed=(sc>0)?0.2*lvl/sc:0;
        if(!d.in[i]){m->SetFixedVariable(o+4+i,Form("bw%d_%d",i+1,b),0);continue;}
        m->SetLimitedVariable(o+4+i,Form("bw%d_%d",i+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*200+1e3); } }
    m->Minimize(); const double* p=m->X();
    std::vector<double> S(NB2,0),D(NB2,0);
    for(int b=0;b<NB2;b++){ SF::BinData&d=SF::B[b]; if(!d.h) continue;
      const double* q=p+2+b*nper;
      double w=std::max(0.05,p[0]),kap=std::max(0.0,p[1]);
      for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i);
        if(x<3.5||x>9.0) continue;
        double bw=0; for(int j=0;j<6;j++) if(d.in[j]) bw+=SF::ampOf(b,j,p)*d.T[j].Eval(x);
        double thr=0.5*(1+std::erf((x-SF::turnOn)/(M_SQRT2*w)));
        double base=d.floor+std::max(0.0,q[2])*(d.gc?std::max(0.0,d.gc->Eval(x)):0.0)
                   +std::max(0.0,q[3])*thr*std::exp(-kap*x)
                   +q[0]*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                   +q[1]*shapeGT(x,d.gsMu+E52d,d.sgL,d.sgR,d.alL,d.alR);
        S[b]+=bw; D[b]+=d.h->GetBinContent(i)-base; } }
    for(int b=0;b<NB2;b++){ if(!SF::B[b].h) continue;
      printf("   %2d  %2.0f-%2.0f %6.0f |   %10.0f    |    %14.0f      | %5.2f\n",
             b,TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,SF::B[b].N,S[b],D[b],D[b]>0?S[b]/D[b]:0); }
    printf("   bin-to-bin ratio of the SUM: ");
    for(int b=0;b+1<NB2;b++) if(S[b]>0&&S[b+1]>0) printf(" %.2f",S[b]/S[b+1]);
    printf("\n");
    delete m; }
}
