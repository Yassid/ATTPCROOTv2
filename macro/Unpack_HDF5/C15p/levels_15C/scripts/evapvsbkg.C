#include "evap.C"
#include "binbuild.C"
// Direct comparison in the SAME theta_CM bins the fit uses: what Ex shape does evaporation
// predict where the background actually lives, and how does it compare with the fitted
// background  B * exp(-kappa * Ex)  (kappa = 0.29 ATTPCROOT, 0.24 Spyral)?
void evapvsbkg(){
  gRandom->SetSeed(20261116);
  setCM(10,5,60);
  double Kb=191.25;
  double Es=EV::Estar(Kb), T=EV::temperature(Es), B=EV::coulombBarrier();
  double pb=std::sqrt(Kb*Kb+2*Kb*EV::M15C);
  TLorentzVector W(0,0,pb,Kb+EV::M15C+EV::MD); TVector3 bv=W.BoostVector();
  Pars A; A.ebeam=Kb; A.mBeam=15.0105993; A.mTarg=2.014102; A.mEject=2.014102; A.mResid=15.0105993;
  A.kcOn=false; A.ebzOn=false;
  const int NCMB=10;
  std::vector<TH1D*> hx(NCMB); std::vector<long> nb(NCMB,0);
  for(int b=0;b<NCMB;b++){ hx[b]=new TH1D(Form("ex%d",b),"",120,-2,22); hx[b]->SetDirectory(nullptr); }
  const int N=20000000; long nacc=0;
  for(int i=0;i<N;i++){
    double eps; while(true){ eps=B+gRandom->Exp(T)+gRandom->Exp(T); if(eps<Es) break; }
    double p=std::sqrt(eps*eps+2*eps*EV::MD);
    double c=-1+2*gRandom->Rndm(), st=std::sqrt(std::max(0.0,1-c*c)), ph=2*M_PI*gRandom->Rndm();
    TLorentzVector d(p*st*std::cos(ph),p*st*std::sin(ph),p*c,eps+EV::MD);
    d.Boost(bv);
    double th=d.Theta()*180/M_PI, ke=d.E()-EV::MD;
    if(ke<0.5||ke>18||th<14||th>86) continue;
    double e=exOf(A,Kb,th,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,Kb,th,ke);
    int b=(int)((tcm-TCM0)/TCMSTEP); if(b<0||b>=NCMB) continue;
    hx[b]->Fill(e); nb[b]++; nacc++;
  }
  printf("  %d evaporated deuterons, %ld (%.3f %%) land in theta_CM 10-60 with the analysis cuts\n",
         N,nacc,100.0*nacc/N);
  printf("\n  bin  theta_CM   evap counts   share   | the background lives in bins 0-1\n");
  for(int b=0;b<NCMB;b++)
    printf("   %2d   %2.0f-%2.0f    %9ld   %6.2f %%\n",b,TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,
           nb[b],nacc?100.0*nb[b]/nacc:0);
  // the Ex shape in bins 0 and 1, against exp(-kappa Ex)
  printf("\n  Ex shape where the background lives (normalised at Ex = 4):\n");
  printf("   Ex          :"); for(double X=3;X<=9.01;X+=1) printf(" %7.1f",X); printf("\n");
  for(int b=0;b<2;b++){
    double n4=hx[b]->GetBinContent(hx[b]->FindBin(4.0)); if(n4<=0) n4=1;
    printf("   evap bin %d  :",b);
    for(double X=3;X<=9.01;X+=1) printf(" %7.3f",hx[b]->GetBinContent(hx[b]->FindBin(X))/n4);
    printf("\n"); }
  for(double kap : {0.29,0.24}){
    printf("   exp(-%.2f Ex):",kap);
    for(double X=3;X<=9.01;X+=1) printf(" %7.3f",std::exp(-kap*(X-4.0)));
    printf("\n"); }
}
