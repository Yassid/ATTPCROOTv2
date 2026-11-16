#include "evap.C"
#include "binbuild.C"
void evaptest(){
  double Kb=191.25;
  double Es=EV::Estar(Kb), T=EV::temperature(Es), B=EV::coulombBarrier();
  printf("  15C + d -> 17N*\n");
  printf("    fusion Q          = %.2f MeV\n",EV::M15C+EV::MD-EV::M17N);
  printf("    E* (compound)     = %.2f MeV\n",Es);
  printf("    T = sqrt(E*/a), a = A/8 = %.3f  ->  T = %.2f MeV\n",17.0/8.0,T);
  printf("    Coulomb barrier for d  = %.2f MeV\n",B);
  // the compound recoils with the beam momentum: velocity and the resulting lab focusing
  double pb=std::sqrt(Kb*Kb+2*Kb*EV::M15C);
  double beta=pb/(Kb+EV::M15C+EV::MD);
  printf("    compound beta     = %.4f   (deuteron at eps = T has beta_cm = %.4f)\n",
         beta,std::sqrt(2*T/ (EV::MD/1.0))*1.0);
  double betaD=std::sqrt(1.0-std::pow(EV::MD/(T+EV::MD),2));
  printf("    deuteron at eps = T:  beta = %.4f  vs compound %.4f  -> %s\n",betaD,beta,
         betaD>beta?"can reach backward lab angles":"FOCUSED FORWARD, cannot reach backward angles");
  double thmax = betaD>beta ? 180.0 : std::asin(betaD/beta)*180/M_PI;
  printf("    max lab angle for an evaporated deuteron at eps = T : %.1f deg\n",thmax);
  printf("\n  the background lives at theta_lab 70-85 deg (theta_CM 10-20), KE ~ 1.4 MeV\n");
}
