#include "binbuild.C"
// Rebuild the theta_CM spectra with 5 deg bins out to 60 deg: 10 angular points instead of 5.
void rebin5(){ setCM(10,5,60); VZLO=0; BB_QUIET=false; binbuild(); }
