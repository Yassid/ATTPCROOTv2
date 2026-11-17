#include "resfitC.C"
// Does dropping 5.98 -- adopted for STABILITY -- cost FIT QUALITY where the data are richest?
// Same configuration otherwise; compare chi2/ndf per bin with and without it.
double chi2b(int b,int withR){
  RH::Bin& d=RH::B[b]; if(!d.h) return -1;
  double s=0; int nd=0;
  for(int i=1;i<=d.h->GetNbinsX();i++){
    double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
    if(x<RH::exLo||x>RH::exHi||n<5) continue;
    double Xt=x-d.gsMu;
    double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
            +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
            +RH::LASTBG[b]*RH::bshape(Xt,RH::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0)
            +(d.in[0]?RH::LASTA[b][0]*d.T[0].Eval(x):0.0)
            +(d.in[1]?RH::LASTA[b][0]*RH::LASTR*d.T[1].Eval(x):0.0)
            +(d.in[3]?RH::LASTA[b][3]*d.T[3].Eval(x):0.0);
    if(withR) m+=(d.in[2]?RH::LASTA[b][2]*d.T[2].Eval(x):0.0);
    if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
  return s/std::max(1,nd-(withR?4:3));
}
void q598(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALON; RH::TIE12=true;
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    double q[2][16]; double L[2];
    for(int set=0;set<2;set++){
      bool keep[6]={true,true,set==1,true,false,false};   // set 1 ADDS 5.98 back
      for(int i=0;i<RH::NR;i++) RH::USE[i]=keep[i];
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
      L[set]=rhRun(CFG[ic],NCM,false);
      for(int b=0;b<NCM;b++) q[set][b]=chi2b(b,set); }
    printf("\n==== %s : chi2/ndf with and without 5.98 ====\n",CFG[ic].nm);
    printf("   %-10s %12s %12s %10s\n","theta_CM","D12+6.68","+5.98","change");
    for(int b=0;b<NCM;b++)
      printf("   %3.0f-%-6.0f %12.2f %12.2f %+10.2f%s\n",
             TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,q[0][b],q[1][b],q[1][b]-q[0][b],
             (q[1][b]<q[0][b]-0.3)?"   <- better with 5.98":"");
    printf("   total -2lnL: %.1f without, %.1f with (10 more parameters)\n",L[0],L[1]);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true; BB_QUIET=true; binbuild();
}
