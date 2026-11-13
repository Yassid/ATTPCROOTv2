#include "binacc.h"
#include "binbuild.C"
// Is the Ex 3-9 window fully accessible in every 5 deg theta_CM bin?  If the kinematic
// acceptance collapses over part of it, the fit is modelling a region the data cannot populate.
void binacc5(){
  setCM(10,5,60);
  for(auto&k:CFG){
    Pars A=parsOf(k); Win w{14,86,k.keLo,k.keHi};
    printf("\n== %s   acceptance over the fit window\n",k.nm);
    printf("   bin  theta_CM |"); for(double X=3;X<=9.01;X+=1) printf("  Ex%4.1f",X); printf("\n");
    for(int b=0;b<10;b++){ double lo=TCM0+b*TCMSTEP, hi=lo+TCMSTEP;
      printf("   %2d   %2.0f-%2.0f   |",b,lo,hi);
      for(double X=3;X<=9.01;X+=1) printf("  %6.3f",accBin(A,k.eb,X,lo,hi,w));
      printf("\n"); }
  }
}
