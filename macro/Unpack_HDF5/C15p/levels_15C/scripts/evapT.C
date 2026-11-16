#include "evap.C"
#include "binbuild.C"
// Does ANY temperature make evaporation look like the background?  The rejection must not depend
// on T = 4.21 MeV, so scan it.  Two things have to come out right simultaneously:
//   (a) the yield must CONCENTRATE in theta_CM bins 0-1 (where the background is);
//   (b) the Ex spectrum there must FALL (the fit gives exp(-0.24..0.29 Ex)).
void evapT(){
  gRandom->SetSeed(7); setCM(10,5,60);
  double Kb=191.25, Es=EV::Estar(Kb), B=EV::coulombBarrier();
  double pb=std::sqrt(Kb*Kb+2*Kb*EV::M15C);
  TLorentzVector W(0,0,pb,Kb+EV::M15C+EV::MD); TVector3 bv=W.BoostVector();
  Pars A; A.ebeam=Kb; A.mBeam=15.0105993; A.mTarg=2.014102; A.mEject=2.014102; A.mResid=15.0105993;
  A.kcOn=false; A.ebzOn=false;
  printf("  T(MeV)   share in theta_CM 10-20   Ex slope there (Ex 3-9)   | data: 100%% and -0.24..-0.29\n");
  for(double T : {1.0,2.0,3.0,4.21,6.0,8.0,12.0,20.0}){
    const int N=3000000; long nlow=0,ntot=0;
    TH1D hx("hx","",60,-2,22); hx.Reset();
    for(int i=0;i<N;i++){
      double eps; int guard=0;
      while(true){ eps=B+gRandom->Exp(T)+gRandom->Exp(T); if(eps<Es||++guard>50) break; }
      if(eps>=Es) continue;
      double p=std::sqrt(eps*eps+2*eps*EV::MD);
      double c=-1+2*gRandom->Rndm(), st=std::sqrt(std::max(0.0,1-c*c)), ph=2*M_PI*gRandom->Rndm();
      TLorentzVector d(p*st*std::cos(ph),p*st*std::sin(ph),p*c,eps+EV::MD);
      d.Boost(bv);
      double th=d.Theta()*180/M_PI, ke=d.E()-EV::MD;
      if(ke<0.5||ke>18||th<14||th>86) continue;
      double e=exOf(A,Kb,th,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,Kb,th,ke);
      if(tcm<10||tcm>60) continue;
      ntot++;
      if(tcm<20){ nlow++; hx.Fill(e); } }
    // log-slope of the Ex spectrum over 3-9 in the 10-20 deg bins
    double sx=0,sy=0,sxx=0,sxy=0; int n=0;
    for(int i=1;i<=hx.GetNbinsX();i++){ double x=hx.GetBinCenter(i),y=hx.GetBinContent(i);
      if(x<3||x>9||y<5) continue; double ly=std::log(y);
      sx+=x; sy+=ly; sxx+=x*x; sxy+=x*ly; n++; }
    double slope = n>2 ? (n*sxy-sx*sy)/(n*sxx-sx*sx) : 0;
    printf("  %5.1f        %7.2f %%                %+8.3f /MeV %s\n",
           T, ntot?100.0*nlow/ntot:0, slope, slope<-0.1?"  <- falls":"  (rises)");
  }
}
