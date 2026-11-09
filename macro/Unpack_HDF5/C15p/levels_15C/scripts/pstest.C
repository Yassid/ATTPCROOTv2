#include "psmodel.h"
void pstest(){
  printf("Sn(15C) from the masses = %.4f MeV   (literature 1.218)\n", PS::SN);
  double Kb=191.25;
  printf("Ex kinematic limit at Ebeam %.2f = %.2f MeV\n", Kb, std::sqrt(PS::sOf(Kb))-PS::MD-PS::M15C);
  printf("\n   Ex      PS density (arb)\n");
  double mx=0; for(double e=0;e<=40;e+=0.1) mx=std::max(mx,PS::dens(e,Kb));
  for(double e=0;e<=24.01;e+=2.0)
    printf("  %5.1f   %.4f\n", e, mx>0? PS::dens(e,Kb)/mx : 0.0);
}
