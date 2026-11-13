#include "psmodel.h"
// Both codes' total angular distributions fall steeply to a minimum near theta_CM 45 deg and
// rise again.  If that is a diffraction minimum, its position fixes qR, hence a radius for an
// assumed transferred angular momentum L (first zero of j_L at qR).
void diffmin(){
  double Kb=191.25;
  double Eb=Kb+PS::M15C, pb=std::sqrt(Kb*Kb+2*Kb*PS::M15C);
  TLorentzVector W(0,0,pb,Eb+PS::MD);
  double s=W.M2(), rs=std::sqrt(s);
  double pcm=PS::omega(s,PS::MD*PS::MD,PS::M15C*PS::M15C)/(2*rs);   // elastic entrance channel
  printf("  CM momentum  k = %.1f MeV/c  =  %.3f fm^-1\n",pcm,pcm/197.327);
  double k=pcm/197.327;
  printf("\n  theta_cm    q = 2k sin(th/2)   (fm^-1)\n");
  for(double th : {15.,25.,35.,45.,55.}){
    double q=2*k*std::sin(th*M_PI/360.0);
    printf("    %4.0f        %8.3f\n",th,q); }
  // first zero of the spherical Bessel j_L
  double z[5]={4.4934,5.7635,6.9879,8.1826,9.3558};   // j_1..j_5 first zeros (j_0 is pi)
  double q45=2*k*std::sin(45*M_PI/360.0);
  printf("\n  if the 45 deg dip is the FIRST diffraction minimum, q = %.3f fm^-1 and\n",q45);
  printf("    L : radius R = (first zero of j_L)/q\n");
  printf("    0 : %6.2f fm\n",M_PI/q45);
  for(int L=1;L<=5;L++) printf("    %d : %6.2f fm\n",L,z[L-1]/q45);
  printf("\n  for reference the channel radius used in the BW penetrability is 1.4*15^(1/3) = %.2f fm,\n",
         1.4*std::pow(15.0,1.0/3.0));
  printf("  and a typical 15C interaction radius is 4-5 fm.\n");
}
