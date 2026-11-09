#include "psmodel.h"
#include "exlib.h"
// Two competing 3-body hypotheses for the continuum, both reconstructed the way the ANALYSIS does
// (two-body missing mass, deuteron mass hypothesis, same cuts):
//   (A) 15C + d -> d + 14C + n     target excitation above Sn, deuteron detected
//   (B) 15C + d -> p + n + 15C     DEUTERON BREAKUP; the PROTON is detected and mis-assigned a
//       deuteron mass, which is what a breakup proton leaking through the gate looks like.
TH1D* gen3(double Kb,const Pars&p,double mA,double mB,double mC,int which,const char* nm,int N=4000000){
  TLorentzVector beam(0,0,std::sqrt(Kb*Kb+2*Kb*PS::M15C),Kb+PS::M15C);
  TLorentzVector W=beam+TLorentzVector(0,0,0,PS::MD);
  double m3[3]={mA,mB,mC}; TGenPhaseSpace g;
  if(!g.SetDecay(W,3,m3)){ printf("%s: SetDecay failed\n",nm); return nullptr; }
  auto* h=new TH1D(nm,";E_{x} [MeV];a.u.",280,-5,23); h->SetDirectory(nullptr);
  long acc=0;
  for(int i=0;i<N;i++){ double w=g.Generate();
    TLorentzVector* d=g.GetDecay(which);
    double ke=d->E()-m3[which], th=d->Theta()*180.0/M_PI;
    if(ke<p.keLo||ke>p.keHi) continue;
    double thc=thCorr(p,th,ke); if(thc<p.thLo||thc>p.thHi) continue;
    double e=exOf(p,Kb,thc,ke);                 // always the DEUTERON hypothesis (exlib Pars)
    if(!std::isfinite(e)) continue; h->Fill(e,w); acc++; }
  printf("  %-26s accepted %ld / %d\n", nm, acc, N);
  return h;
}
void psalt(){
  Pars A; A.ebeam=191.25; A.keLo=0.5; A.keHi=18; A.thLo=0; A.thHi=180;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;
  const double MP=1.00727646*PS::U;
  TH1D* hA=gen3(191.25,A,PS::MD,PS::M14C,PS::MN,0,"A_d_14C_n");        // detect the deuteron
  TH1D* hB=gen3(191.25,A,MP,PS::MN,PS::M15C,0,"B_breakup_proton");     // detect the proton
  TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* hd=(TH1D*)f.Get("ATTPCROOT"); hd->SetDirectory(nullptr);
  auto nrm=[&](TH1D* h){ int b=h->FindBin(10.0); double d=hd->GetBinContent(hd->FindBin(10.0));
    if(h->GetBinContent(b)>0) h->Scale(d/h->GetBinContent(b)); };
  if(hA) nrm(hA); if(hB) nrm(hB);
  printf("\n   %6s %9s %11s %11s\n","Ex","data","(A) d+14C+n","(B) breakup p");
  for(double e=2;e<=20.01;e+=2.0){ int b=hd->FindBin(e);
    printf("   %6.1f %9.0f %11.0f %11.0f\n", e, hd->GetBinContent(b),
      hA?hA->GetBinContent(hA->FindBin(e)):0, hB?hB->GetBinContent(hB->FindBin(e)):0); }
  TFile o("/home/yassid/Data/C15d_statefit/ps_alt.root","RECREATE");
  if(hA)hA->Write("A"); if(hB)hB->Write("B"); o.Close();
}
