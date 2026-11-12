#include "binacc.h"
#include "binbuild.C"
void kedbg(){
  Pars A=parsOf(CFG[0]);
  printf("labOf: theta_lab and KE vs Ex, at fixed theta_CM\n");
  for(double cm : {15.0,25.0,35.0,45.0}){
    printf("  theta_CM %4.0f :",cm);
    for(double Ex=0;Ex<=20.01;Ex+=2.5){ double thl,kk; labOf(A,191.25,Ex,cm,thl,kk);
      printf("  Ex%4.1f(%5.1fdeg,%5.1fMeV)",Ex,thl,kk); }
    printf("\n"); }
}
