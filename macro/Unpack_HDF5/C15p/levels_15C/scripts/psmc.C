#include "psmodel.h"
#include "exlib.h"
// Phase-space Monte Carlo WITH the analysis acceptance.
//   generate 15C + d -> d + 14C + n uniformly in 3-body phase space (TGenPhaseSpace),
//   take the DEUTERON's lab KE and theta, apply the SAME cuts the data passed,
//   reconstruct Ex the SAME way (two-body missing mass from the deuteron alone),
//   histogram -> the predicted continuum shape.
// The reconstructed Ex of a 3-body event is NOT the generated invariant mass once the cuts bite,
// which is exactly why this has to be simulated rather than taken from the analytic density.
TH1D* psmc(double Kb, const Pars& p, const char* nm, int N=4000000, unsigned seed=12345){
  TRandom3 rng(seed);
  TLorentzVector beam(0,0,std::sqrt(Kb*Kb+2*Kb*PS::M15C), Kb+PS::M15C);
  TLorentzVector targ(0,0,0,PS::MD);
  TLorentzVector W = beam + targ;
  double masses[3] = {PS::MD, PS::M14C, PS::MN};
  TGenPhaseSpace gen;
  if(!gen.SetDecay(W,3,masses)){ printf("SetDecay failed\n"); return nullptr; }
  auto* h=new TH1D(nm,";E_{x} [MeV];a.u.",280,-5,23); h->SetDirectory(nullptr);
  long acc=0;
  for(int i=0;i<N;i++){
    double w=gen.Generate();                       // phase-space weight
    TLorentzVector* d=gen.GetDecay(0);             // the detected deuteron
    double ke=d->E()-PS::MD, th=d->Theta()*180.0/M_PI;
    if(ke<p.keLo||ke>p.keHi) continue;
    double thc=thCorr(p,th,ke);
    if(thc<p.thLo||thc>p.thHi) continue;
    // reconstruct Ex exactly as the analysis does, from the deuteron alone
    double e=exOf(p,Kb,thc,ke);
    if(!std::isfinite(e)) continue;
    h->Fill(e,w); acc++; }
  printf("  %-12s generated %d, accepted %ld (%.2f %%)\n", nm, N, acc, 100.0*acc/N);
  return h;
}
