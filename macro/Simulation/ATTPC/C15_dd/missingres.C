#include "/home/yassid/Data/C15d_statefit/excm.h"
// The simulation reproduces the efficiency but is ~2x too GOOD in resolution at Ex = 0.
// How much width is missing, and could a beam-energy spread account for it?
void missingres(){
  // measured (lineshape.C on dd_ex.root) vs simulated (response_C15dd.C) at Ex = 0
  struct R{const char* nm; double fL,fR;};
  R data[2]={{"ATTPCROOT",0.752,1.142},{"Spyral",0.497,0.523}};
  double sL=0.234,sR=0.279;                       // simulation
  double sigSim=0.5*(sL+sR)/2.355;
  printf("  simulated sigma at Ex=0 : %.3f MeV  (FWHM %.3f / %.3f)\n\n",sigSim,sL,sR);
  // dEx / dEbeam, numerically, at the g.s. for a typical accepted deuteron
  Pars P; P.ebeam=191.25; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=false;
  const double md=2.0141017781*931.49401;
  // take theta_lab 75 deg, the bulk of the data; KE from elastic kinematics at that angle
  double th=75.0;
  double ke0=0; { // solve Ex=0 for KE
    double lo=0.2,hi=40; for(int i=0;i<60;i++){ double m=0.5*(lo+hi);
      double e=exOf(P,191.25,th,m); if(!std::isfinite(e)||e>0) lo=m; else hi=m; } ke0=0.5*(lo+hi); }
  double d=0.5, e1=exOf(P,191.25-d,th,ke0), e2=exOf(P,191.25+d,th,ke0);
  double dExdEb=(e2-e1)/(2*d);
  printf("  at theta_lab %.0f deg, KE %.2f MeV :  dEx/dEbeam = %+.4f\n",th,ke0,dExdEb);
  printf("\n  %-10s  measured sigma   missing (quadrature)   implied beam spread sigma(E)\n","code");
  for(auto&r:data){ double sg=0.5*(r.fL+r.fR)/2.355;
    double miss=std::sqrt(std::max(0.0,sg*sg-sigSim*sigSim));
    printf("  %-10s   %.3f MeV        %.3f MeV              %.2f MeV  (%.2f %% of 191)\n",
           r.nm,sg,miss,std::fabs(dExdEb)>1e-6?miss/std::fabs(dExdEb):0,
           std::fabs(dExdEb)>1e-6?100*miss/std::fabs(dExdEb)/191.25:0); }
}
