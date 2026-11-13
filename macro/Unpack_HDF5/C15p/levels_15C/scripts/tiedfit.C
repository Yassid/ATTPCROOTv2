#include "simfit.C"
// Fit RESOLVABLE STRUCTURES: the two unresolvable doublets tied with a shared internal ratio.
//   E1 = BW1 + r1*BW2 (4.43+4.93)   E2 = BW3 (5.98)   E3 = BW4 + r2*BW5 (6.68+7.38)   E4 = BW6
// 4 amplitudes per bin + 2 shared ratios = 22 free, against 30 for six free amplitudes.
double tiedRun(const Cfg& k,bool talk){
  SF::B.clear(); SF::UCM.clear();
  for(int ib=0;ib<SF::NBIN;ib++) sfPrep(k,ib);
  for(int ib=0;ib<SF::NBIN;ib++) SF::UCM.push_back(std::cos((TCM0+(ib+0.5)*TCMSTEP)*M_PI/180.0));
  SF::SMOOTH=false; SF::TIED=true;
  int npar=2+SF::NBIN*8+2;
  ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
  m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000);
  m->SetTolerance(0.1); m->SetPrintLevel(0); m->SetStrategy(1);
  ROOT::Math::Functor fnc(&SF::fcn,npar); m->SetFunction(fnc);
  m->SetLimitedVariable(0,"bkg_w",1.0,0.02,0.05,3.0);
  m->SetLimitedVariable(1,"bkg_kappa",0.3,0.01,0.0,2.0);
  for(int b=0;b<SF::NBIN;b++){ SF::BinData&d=SF::B[b]; if(!d.h) continue;
    double pk=d.h->GetBinContent(d.h->GetMaximumBin()); double lvl=0;int nl=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){double x=d.h->GetBinCenter(i);
      if(x>6&&x<8){lvl+=d.h->GetBinContent(i);nl++;}} lvl=nl?lvl/nl:1;
    int o=2+b*8;
    m->SetLimitedVariable(o+0,Form("gs_%d",b),pk,pk*1e-3,0,pk*3);
    m->SetLimitedVariable(o+1,Form("e52_%d",b),pk*0.05,pk*1e-4,0,pk);
    m->SetLimitedVariable(o+2,Form("cont_%d",b),1e-3,1e-5,0,1e3);
    m->SetLimitedVariable(o+3,Form("bkgB_%d",b),0.5*lvl,lvl*1e-3,0,lvl*20);
    for(int e=0;e<4;e++){ int a=SF::EMAP[e][0];
      double sc=d.in[a]?d.T[a].Eval((EPP[a]-PPGS)+d.gsMu):0;
      double seed=(sc>0)?0.2*lvl/sc:0;
      m->SetLimitedVariable(o+4+e,Form("E%d_%d",e+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*200+1e3); } }
  m->SetLimitedVariable(2+SF::NBIN*8+0,"r1_4.93/4.43",1.0,0.02,0.0,10.0);
  m->SetLimitedVariable(2+SF::NBIN*8+1,"r2_7.38/6.68",0.4,0.02,0.0,10.0);
  m->Minimize(); const double* p=m->X();
  if(talk){
    printf("\n================ %s   TIED doublets (%d parameters) ================\n",k.nm,npar);
    printf("   shared: bkg w %.3f  kappa %.4f | r1(4.93/4.43) = %.3f   r2(7.38/6.68) = %.3f   -2lnL %.1f\n",
           p[0],p[1],p[2+SF::NBIN*8+0],p[2+SF::NBIN*8+1],m->MinValue());
    printf("   bin  cm      N   chi2/ndf |   E1(4.4+4.9)   E2(5.98)   E3(6.7+7.4)   E4(8.28)\n");
    std::vector<double> A[4];
    double tot=0; int ntot=0;
    for(int b=0;b<SF::NBIN;b++){ SF::BinData&d=SF::B[b]; if(!d.h) continue;
      const double* q=p+2+b*8;
      double c2=0;int nd=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
        if(x<SF::exLo||x>SF::exHi||n<10) continue;
        double mm=SF::model(b,x,p); c2+=(n-mm)*(n-mm)/n; nd++; }
      tot+=c2; ntot+=nd;
      printf("   %2d  %2.0f-%2.0f %6.0f  %6.2f  |",b,TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,d.N,c2/std::max(1,nd-8));
      for(int e=0;e<4;e++){ int a=SF::EMAP[e][0];
        double ph=d.in[a]?q[4+e]*d.T[a].Eval((EPP[a]-PPGS)+d.gsMu):0;
        A[e].push_back(ph); printf(" %12.1f",ph); }
      printf("\n"); }
    printf("   overall chi2/ndf = %.3f\n",tot/std::max(1,ntot-npar));
    printf("   bin-to-bin |ratio| :");
    for(int e=0;e<4;e++){ double worst=0;
      for(size_t b=0;b+1<A[e].size();b++){ double x=A[e][b],y=A[e][b+1];
        if(x>1&&y>1) worst=std::max(worst,x>y?x/y:y/x); else if(x>5||y>5) worst=99; }
      printf("   E%d %5.1f",e+1,worst); }
    printf("    (free-amplitude version gave 5.1, 5.2, 99, 3.2, 99, 99)\n");
  }
  double v=m->MinValue(); SF::TIED=false; delete m; return v;
}
void tiedfit(){ for(auto&k:CFG) tiedRun(k,true); }
