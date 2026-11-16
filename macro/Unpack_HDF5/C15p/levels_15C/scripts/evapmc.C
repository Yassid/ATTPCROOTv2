#include "evap.C"
#include "binbuild.C"
// Full evaporation Monte Carlo: sample eps from the surface Maxwellian above the barrier,
// emit isotropically in the compound frame, boost to the lab, apply the analysis cuts and
// reconstruct Ex with the SAME two-body deuteron hypothesis the data uses.
void evapmc(){
  gRandom->SetSeed(20261116);
  double Kb=191.25;
  double Es=EV::Estar(Kb), T=EV::temperature(Es), B=EV::coulombBarrier();
  double pb=std::sqrt(Kb*Kb+2*Kb*EV::M15C);
  TLorentzVector W(0,0,pb,Kb+EV::M15C+EV::MD);           // the compound 4-vector
  TVector3 bv=W.BoostVector();
  Pars A; A.ebeam=Kb; A.mBeam=15.0105993; A.mTarg=2.014102; A.mEject=2.014102; A.mResid=15.0105993;
  A.kcOn=false; A.ebzOn=false;
  const int N=4000000;
  TH1D hAll("hAll",";theta_lab",90,0,90), hAcc("hAcc",";theta_lab",90,0,90);
  TH1D hEx("hEx",";Ex",120,-2,22), hKE("hKE",";KE",80,0,40);
  TH2D hM("hM",";theta_lab;KE",90,0,90,80,0,40);
  long nAcc=0, nBack=0;
  for(int i=0;i<N;i++){
    // surface Maxwellian above the barrier, by rejection
    double eps;
    while(true){ eps=B+gRandom->Exp(T)+gRandom->Exp(T);   // (eps-B) exp(-(eps-B)/T) = Gamma(2,T)
                 if(eps<Es) break; }
    double p=std::sqrt(eps*eps+2*eps*EV::MD);
    double c=-1+2*gRandom->Rndm(), st=std::sqrt(std::max(0.0,1-c*c)), ph=2*M_PI*gRandom->Rndm();
    TLorentzVector d(p*st*std::cos(ph),p*st*std::sin(ph),p*c,eps+EV::MD);
    d.Boost(bv);
    double th=d.Theta()*180/M_PI, ke=d.E()-EV::MD;
    hAll.Fill(th);
    if(th>65) nBack++;
    // the analysis cuts
    if(ke<0.5||ke>18) continue;
    if(th<14||th>86) continue;
    hAcc.Fill(th); hKE.Fill(ke); hM.Fill(th,ke); nAcc++;
    double e=exOf(A,Kb,th,ke);
    if(std::isfinite(e)) hEx.Fill(e);
  }
  printf("  E* = %.1f MeV, T = %.2f MeV, barrier = %.2f MeV,  %d evaporated deuterons\n",Es,T,B,N);
  printf("  reaching theta_lab > 65 deg : %ld  (%.3f %%)\n",nBack,100.0*nBack/N);
  printf("  passing the analysis cuts   : %ld  (%.2f %%)\n",nAcc,100.0*nAcc/N);
  printf("\n  theta_lab distribution of ALL evaporated deuterons:\n   ");
  for(int b=0;b<9;b++){ double s=hAll.Integral(hAll.FindBin(b*10+0.5),hAll.FindBin(b*10+9.5));
    printf(" %d-%d:%5.2f%%",b*10,b*10+10,100.0*s/N); }
  printf("\n  DATA for comparison (deuteron gate):\n");
  printf("    0-10:0.0%% 10-20:1.5%% 20-30:2.8%% 30-40:3.4%% 40-50:4.3%% 50-60:6.5%% 60-70:7.2%% 70-80:35.4%% 80-90:38.8%%\n");
  printf("\n  reconstructed Ex of the accepted evaporated deuterons:\n   ");
  for(double X=0;X<=20.01;X+=2) printf(" %4.0f:%5.0f",X,hEx.GetBinContent(hEx.FindBin(X)));
  printf("\n");
  TFile fo("/home/yassid/Data/C15d_statefit/evap_mc.root","RECREATE");
  hAll.Write(); hAcc.Write(); hEx.Write(); hKE.Write(); hM.Write(); fo.Close();
}
