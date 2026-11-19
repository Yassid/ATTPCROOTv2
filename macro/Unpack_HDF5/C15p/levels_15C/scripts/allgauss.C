#include "two1015.C"
// Every angular bin, both codes, both background cases.  Peaks are PURE GAUSSIANS -- no tails.
// chi2 is reported split: the elastic peak region (|Ex| < 2.5) and the RESONANCE region above it.
// They are effectively independent, since the g.s. contributes 0.00-0.20 counts at Ex = 4.
void allgauss(){ for(int b=0;b<NCM;b++) two1015(b); }
