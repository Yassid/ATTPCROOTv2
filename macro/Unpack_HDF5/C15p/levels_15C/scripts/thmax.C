#include "psmodel.h"
// In inverse kinematics the scattered deuteron has a MAXIMUM lab angle whenever its CM speed
// drops below the CM speed itself.  At Ex = 0 the two are equal and theta_max = 90 deg (the
// classic elastic-recoil result).  As Ex grows, beta_d(cm) falls and theta_max closes in.
// The detector only sees theta >~ 65 deg (below that the recoil KE exceeds the 18 MeV cut),
// so there is a HARD kinematic ceiling on the Ex reachable in this channel.
void thmax(){
  double Kb=191.25;
  double Eb=Kb+PS::M15C, pb=std::sqrt(Kb*Kb+2*Kb*PS::M15C);
  TLorentzVector W(0,0,pb,Eb+PS::MD);
  double s=W.M2(), rs=std::sqrt(s);
  double bcm=W.BoostVector().Mag();
  printf("beam 15C %.2f MeV on d:  sqrt(s) = %.1f   beta_cm = %.5f\n",Kb,rs,bcm);
  printf("kinematic limit  Ex_max = sqrt(s) - M(15C) - M(d) = %.2f MeV\n\n", rs-PS::M15C-PS::MD);
  printf("   Ex    beta_d(cm)   theta_max   KE at 70deg   KE at 80deg   reachable at theta>70?\n");
  for(double X=0;X<=20.01;X+=1.0){
    double M4=PS::M15C+X, M3=PS::MD;
    if(rs<=M3+M4){ printf("  %5.1f    (above the kinematic limit)\n",X); continue; }
    double pcm=PS::omega(s,M3*M3,M4*M4)/(2*rs), E3=std::sqrt(pcm*pcm+M3*M3);
    double bd=pcm/E3;
    double thmax = (bd>=bcm)? 90.0 : std::asin(bd/bcm)*180.0/M_PI;
    // KE at a given lab angle: scan cos(theta_cm)
    auto keAt=[&](double thl)->double{ double best=-1, dmin=1e9;
      for(int i=0;i<4000;i++){ double c=-1.0+2.0*(i+0.5)/4000, st=std::sqrt(std::max(0.0,1-c*c));
        TLorentzVector d(pcm*st,0,pcm*c,E3); d.Boost(W.BoostVector());
        double t=d.Theta()*180.0/M_PI; if(std::fabs(t-thl)<dmin){dmin=std::fabs(t-thl); best=d.E()-M3;} }
      return dmin<1.0? best : -1; };
    double k70=keAt(70.0), k80=keAt(80.0);
    printf("  %5.1f     %.5f     %6.2f  %12s %13s   %s\n",X,bd,thmax,
           k70>0?Form("%.2f",k70):"--", k80>0?Form("%.2f",k80):"--",
           thmax>70?"yes":"NO -- forbidden");
  }
}
