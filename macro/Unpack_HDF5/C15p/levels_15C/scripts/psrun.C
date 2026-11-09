#include "psmc.C"
void psrun(){
  gStyle->SetOptStat(0);
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.keLo=0.5; A.keHi=18; A.thLo=0; A.thHi=180;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;   // ebz needs vz, which the MC has no notion of
  Pars S=A; S.ebeam=195; S.keLo=0; S.keHi=20; S.kcDenom=2535.2; S.kcPivot=14.5;
  printf("3-body phase space + acceptance:\n");
  TH1D* ha=psmc(191.25,A,"PS_ATTPCROOT");
  TH1D* hs=psmc(195.00,S,"PS_Spyral");
  TFile o("/home/yassid/Data/C15d_statefit/ps_mc.root","RECREATE");
  if(ha){ha->Write("PS_ATTPCROOT");} if(hs){hs->Write("PS_Spyral");}
  o.Close();
  printf("\n   Ex    PS+acc (ATTPCROOT, normalised to its max)\n");
  double mx=ha?ha->GetMaximum():1;
  for(double e=0;e<=22.01;e+=2.0){ int b=ha->FindBin(e);
    printf("  %5.1f   %.3f\n", e, mx>0? ha->GetBinContent(b)/mx : 0.0); }
}
