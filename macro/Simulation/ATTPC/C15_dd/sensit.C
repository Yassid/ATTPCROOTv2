#include "/home/yassid/Data/C15d_statefit/excm.h"
// Ex is insensitive to Ebeam (dEx/dEbeam = 0.0155), so a beam-energy spread cannot supply the
// missing resolution -- it would need 6-13 % of 191 MeV.  Which MEASURED quantity does, then?
void sensit(){
  Pars P; P.ebeam=191.25; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=false;
  printf("  th_lab   KE     dEx/dtheta   dEx/dKE   | sigma(theta) needed for the missing width\n");
  printf("   (deg)  (MeV)   (MeV/deg)   (MeV/MeV) |  ATTPCROOT 0.387     Spyral 0.187  (MeV)\n");
  for(double th : {65.0,70.0,75.0,80.0,85.0}){
    double ke0=0; { double lo=0.2,hi=40;
      for(int i=0;i<60;i++){ double m=0.5*(lo+hi); double e=exOf(P,191.25,th,m);
        if(!std::isfinite(e)||e>0) lo=m; else hi=m; } ke0=0.5*(lo+hi); }
    double dt=0.25, dk=0.05;
    double dEdth=(exOf(P,191.25,th+dt,ke0)-exOf(P,191.25,th-dt,ke0))/(2*dt);
    double dEdke=(exOf(P,191.25,th,ke0+dk)-exOf(P,191.25,th,ke0-dk))/(2*dk);
    printf("   %4.0f   %5.2f   %+9.4f   %+8.4f  |   %8.3f deg      %8.3f deg\n",
           th,ke0,dEdth,dEdke,
           std::fabs(dEdth)>1e-9?0.387/std::fabs(dEdth):0,
           std::fabs(dEdth)>1e-9?0.187/std::fabs(dEdth):0);
  }
  printf("\n  and if it were all KE instead:\n");
  double th=75.0, ke0=5.26, dk=0.05;
  double dEdke=(exOf(P,191.25,th,ke0+dk)-exOf(P,191.25,th,ke0-dk))/(2*dk);
  printf("   sigma(KE) needed:  ATTPCROOT %.3f MeV (%.1f %% of %.2f)   Spyral %.3f MeV (%.1f %%)\n",
         0.387/std::fabs(dEdke),100*0.387/std::fabs(dEdke)/ke0,ke0,
         0.187/std::fabs(dEdke),100*0.187/std::fabs(dEdke)/ke0);
}
