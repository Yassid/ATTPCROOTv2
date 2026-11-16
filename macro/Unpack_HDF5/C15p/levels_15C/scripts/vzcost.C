#include "resfitC.C"
// What does VZLO = 200 actually cost, and does it change the answer?
// The vz 90-180 mm excess is instrumental, so if it is there it has been contaminating the
// RESONANCE region too, not just the continuum -- that is the real reason to cut it.
void vzcost(){
  setCM(10,5,60);
  double cuts[3]={0,100,200};
  printf("\n================ statistics cost of the low-vz cut ================\n");
  for(int ic=0;ic<2;ic++){
    printf("   %-10s ",CFG[ic].nm);
    double N0=0;
    for(size_t j=0;j<3;j++){
      VZLO=cuts[j]; BB_QUIET=true; binbuild();
      TFile f("dd_cmbins.root"); double N=0;
      for(int b=0;b<NCM;b++){ TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",CFG[ic].nm,b));
        if(h) N+=h->GetEntries(); }
      f.Close();
      if(j==0){ N0=N; printf("VZLO=0: %8.0f   ",N); }
      else printf("VZLO=%3.0f: %8.0f (%.1f%% kept)   ",cuts[j],N,100.0*N/N0); }
    printf("\n"); }
  VZLO=0; BB_QUIET=true; binbuild();
}
