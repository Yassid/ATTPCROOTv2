#include "one1015.C"
// ANGULAR DISTRIBUTIONS from the final configuration.
//   peaks: pure Gaussians (centroid + common width fitted)
//   background: p3* q (Ex-Sn)^a exp(-k(Ex-Sn)), free PER ANGLE, spanning g.s. to 10 MeV
//   states fixed; pair ratios from the global fit
// Every bin is fitted TWICE -- VZLO 300 and 350, about a 5 % change in the sample -- and the
// spread enters the error.  The shape parameter `a` is SCANNED (Migrad will not find its global
// minimum), and a bin whose two passes disagree by more than a factor 2 is flagged, not plotted.
static double YL[2][2][16][2];     // [pass][code][bin][P1,P2]
static double QC[2][2][16];        // resonance-region chi2/ndf
static double SH[2][16][2];        // [code][bin][a,k] from pass 0
void angfinal(){
  setCM(10,5,60);
  // ⚠ FIT from Ex = 2.0, not from -1.5.  The background is still DEFINED from Sn -- there is no
  // gap in the model -- but the elastic peak is excluded from the FIT.  A pure Gaussian describes
  // that peak badly (peak chi2 ~ 14-38 forward), and with thousands of counts in it that misfit
  // dragged ATTPCROOT's background steep: a=2.00 k=2.03, chi2R 4.31, and P1 inflated to 55.5
  // against Spyral's 18.4.  Fitting from 2.0 gives the two codes IDENTICAL shape parameters
  // (a=1.50, k=0.64) and brings P1 to 20.6 vs 15.4.
  B1::exLo=2.0;
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int pass=0;pass<2;pass++){
    for(int ic=0;ic<2;ic++){
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
      VZLO = pass? 350 : 300; BB_QUIET=true; binbuild();
      RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
      for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
      B1::EB=CFG[ic].eb;
      for(int b=0;b<NCM;b++){
        RK::Bin& d=RK::B[b]; if(!d.h){ YL[pass][ic][b][0]=YL[pass][ic][b][1]=-1; continue; }
        B1::BIN=b; RK::ONLYBIN=b;
        Pars A=parsOf(CFG[ic]); A.vzLo=VZLO; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
        B1::ACC=new TGraph();
        for(int i=0;i<=460;i++){ double X=i*0.05;
          B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
        const int npar=3+2+4+2; B1::NPAR=npar; B1::NSH=2; B1::FORM=3;
        double pk=d.h->GetBinContent(d.h->GetMaximumBin());
        std::vector<double> best(npar,0); double bestL=1e18;
        const double AG[9]={0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0};
        for(double av:AG){
          auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
          m->SetMaxFunctionCalls(500000); m->SetTolerance(0.01);
          m->SetPrintLevel(0); m->SetStrategy(2);
          ROOT::Math::Functor fn(&B1::fcn,npar); m->SetFunction(fn);
          m->SetLimitedVariable(0,"Ags",pk,pk*0.01,0,pk*10);
          m->SetLimitedVariable(1,"A52",0.1*pk,pk*0.001,0,pk*5);
          m->SetLimitedVariable(2,"B",50,1,0,1e7);
          m->SetFixedVariable(3,"a",av);
          m->SetLimitedVariable(4,"k",0.5,0.02,0.0,9.0);
          for(int g=0;g<4;g++) m->SetLimitedVariable(5+g,Form("A%d",g),20,0.5,0,1e5);
          m->SetLimitedVariable(npar-2,"mu",d.gsMu,0.005,d.gsMu-0.6,d.gsMu+0.6);
          m->SetLimitedVariable(npar-1,"sig",0.4,0.005,0.08,1.5);
          m->Minimize();
          if(m->MinValue()<bestL){ bestL=m->MinValue(); best.assign(m->X(),m->X()+npar); }
          delete m; }
        // ⚠ the scan FIXES a at a grid point, so a was never actually fitted -- every value
        // came back exactly on the grid.  Having found the right basin, release it and refit.
        { auto* mf=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
          mf->SetMaxFunctionCalls(500000); mf->SetTolerance(0.01);
          mf->SetPrintLevel(0); mf->SetStrategy(2);
          ROOT::Math::Functor fnf(&B1::fcn,npar); mf->SetFunction(fnf);
          mf->SetLimitedVariable(0,"Ags",best[0],std::max(1e-3,best[0]*0.01),0,pk*10);
          mf->SetLimitedVariable(1,"A52",best[1],std::max(1e-3,best[1]*0.01+1e-3),0,pk*5);
          mf->SetLimitedVariable(2,"B",best[2],std::max(1e-3,best[2]*0.01),0,1e7);
          // seed just OFF the grid point so a starts interior, never on a limit
          mf->SetLimitedVariable(3,"a",std::max(0.02,best[3]),0.02,0.0,9.0);
          mf->SetLimitedVariable(4,"k",best[4],0.02,0.0,9.0);
          for(int g=0;g<4;g++) mf->SetLimitedVariable(5+g,Form("A%d",g),best[5+g],
                                   std::max(1e-3,best[5+g]*0.01+1e-3),0,1e5);
          mf->SetLimitedVariable(npar-2,"mu",best[npar-2],0.005,d.gsMu-0.6,d.gsMu+0.6);
          mf->SetLimitedVariable(npar-1,"sig",best[npar-1],0.005,0.08,1.5);
          mf->Minimize();
          if(mf->MinValue()<bestL){ bestL=mf->MinValue(); best.assign(mf->X(),mf->X()+npar); }
          delete mf; }
        const double* p=best.data();
        // peak heights of the two tied pairs
        auto grpH=[&](int g){ double mx=0; const double* a2=p+5;
          for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==g&&d.in[r]){
            double mu=(RK::EPP6[r]-RK::PPGS6)+d.gsMu, sm=0;
            for(int q=0;q<RK::NR;q++) if(RK::GRP[q]==g&&d.in[q])
              sm+=((RK::RIDX[q]>=0)?RK::FIXRAT[RK::RIDX[q]]:1.0)*d.T[q].Eval(mu);
            mx=std::max(mx,std::max(0.0,a2[g])*sm); }
          return mx; };
        YL[pass][ic][b][0]=grpH(0); YL[pass][ic][b][1]=grpH(1);
        if(!pass){ SH[ic][b][0]=p[3]; SH[ic][b][1]=p[4]; }
        double s2R=0; int ndR=0;
        for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
          if(x<B1::PKHI||x>B1::exHi) continue;
          double mm=B1::model(x,p); if(mm<1e-9) mm=1e-9;
          s2R+=(n-mm)*(n-mm)/std::max(mm,1.0); ndR++; }
        QC[pass][ic][b]=s2R/std::max(1.0,ndR-7.0);
        RK::ONLYBIN=-1; } } }
  // ---------- report ----------
  const char* nm[2]={"ATTPCROOT","Spyral"};
  for(int ic=0;ic<2;ic++){
    printf("\n================ %s : angular distributions, final configuration ================\n",nm[ic]);
    printf("   %-9s %8s %8s %9s | %8s %8s %9s | %7s %7s %6s\n",
           "theta_CM","P1","refit","stab","P2","refit","stab","chi2R","a","k");
    for(int b=0;b<NCM;b++){
      double p10=YL[0][ic][b][0],p11=YL[1][ic][b][0];
      double p20=YL[0][ic][b][1],p21=YL[1][ic][b][1];
      auto st=[&](double a,double c){ if(a<=0.5||c<=0) return std::string("  --  ");
        double r=c/0.94/a; char s[16]; snprintf(s,16,"%5.2f%s",r,(r<0.5||r>2.0)?"!":" ");
        return std::string(s); };
      printf("   %3.0f-%-5.0f %8.1f %8.1f %9s | %8.1f %8.1f %9s | %7.2f %7.2f %6.2f\n",
             TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,p10,p11,st(p10,p11).c_str(),
             p20,p21,st(p20,p21).c_str(),QC[0][ic][b],SH[ic][b][0],SH[ic][b][1]); }
    printf("   ! = the two passes disagree by more than a factor 2 -- not a measurement\n");
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
