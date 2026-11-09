#include "psmc.C"
// If the continuum IS 3-body phase space, then data/PS is the reconstruction efficiency, and it
// must be a smooth monotonic function of the DEUTERON ENERGY (short, slow tracks reconstruct
// worse). If instead it is ragged or non-monotonic, the phase-space interpretation fails.
void effcheck(){
  gStyle->SetOptStat(0);
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.keLo=0.5; A.keHi=18; A.thLo=0; A.thHi=180;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;
  // mean deuteron KE vs reconstructed Ex, inside the acceptance (from the same MC)
  TRandom3 rng(7);
  TLorentzVector beam(0,0,std::sqrt(191.25*191.25+2*191.25*PS::M15C),191.25+PS::M15C);
  TLorentzVector W=beam+TLorentzVector(0,0,0,PS::MD);
  double m3[3]={PS::MD,PS::M14C,PS::MN}; TGenPhaseSpace gen; gen.SetDecay(W,3,m3);
  auto* hps=new TH1D("hps","",56,-5,23); hps->SetDirectory(nullptr);
  auto* hke=new TProfile("hke","",56,-5,23); hke->SetDirectory(nullptr);
  for(int i=0;i<4000000;i++){ double w=gen.Generate(); TLorentzVector* d=gen.GetDecay(0);
    double ke=d->E()-PS::MD, th=d->Theta()*180.0/M_PI;
    if(ke<A.keLo||ke>A.keHi) continue;
    double e=exOf(A,191.25,thCorr(A,th,ke),ke);
    if(!std::isfinite(e)) continue; hps->Fill(e,w); hke->Fill(e,ke,w); }
  // the data
  TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* h0=(TH1D*)f.Get("ATTPCROOT"); h0->SetDirectory(nullptr);
  auto* hd=new TH1D("hd","",56,-5,23); hd->SetDirectory(nullptr);
  for(int b=1;b<=h0->GetNbinsX();b++) hd->Fill(h0->GetBinCenter(b),h0->GetBinContent(b));
  printf("   %6s %10s %10s %10s %10s\n","Ex","data","PS+acc","data/PS","<KE_d>");
  double norm=0; { int b=hd->FindBin(6.0); norm = hps->GetBinContent(b)>0 ? hd->GetBinContent(b)/hps->GetBinContent(b) : 1; }
  for(double e=2;e<=20.01;e+=1.0){ int b=hd->FindBin(e);
    double D=hd->GetBinContent(b), P=hps->GetBinContent(b)*norm;
    printf("   %6.1f %10.0f %10.0f %10.3f %10.2f\n", e, D, P, P>0?D/P:0, hke->GetBinContent(b)); }
}
