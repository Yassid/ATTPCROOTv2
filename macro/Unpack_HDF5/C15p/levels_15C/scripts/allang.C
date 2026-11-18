#include "two1015.C"
// OVERNIGHT: every angular bin, both codes, both background cases, and BOTH peak shapes.
//   peaks GAUSS  : pure Gaussian g.s. and 5/2+, centroid and a common width fitted
//   peaks GT     : the bifurcated Gaussian + exponential tails used so far
// Everything else as agreed: fit g.s. to 10 MeV with no gap, states fixed, pair ratios from the
// global fit, background either free (2 shape par) or constrained (1).
void allang(){
  for(int g=1;g>=0;g--){            // GAUSS first, then GT for comparison
    B1::GAUSS = (g==1);
    printf("\n\n##################################################################\n");
    printf("###  PEAK SHAPE: %s\n", B1::GAUSS? "PURE GAUSSIAN (g.s. and 5/2+)"
                                             : "bifurcated Gaussian + exponential tails");
    printf("##################################################################\n");
    for(int b=0;b<NCM;b++) two1015(b);
  }
}
