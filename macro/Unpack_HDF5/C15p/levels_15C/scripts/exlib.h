#pragma once
// Reproduces the explorer's physics EXACTLY (page lines 408-575):
//   U = 931.49401 ; M_i = amu_i * U
//   theta_corr = theta - (360/kcDenom)*(KE + keOff - kcPivot)      [deg]
//   Ebeam(vz)  = ebeam - (ebzSlope/1000)*(vz - ebzPivot)           [MeV, vz mm]
//   Ex from the same kine2b expression, so numbers here match what is on screen.
#include <cmath>
static const double U_AMU = 931.49401;
struct Pars {
  double mBeam=15.010599, mTarg=2.014102, mEject=2.014102, mResid=15.010599;
  double ebeam=198, chi2=5, thLo=0, thHi=180, keLo=0, keHi=1000;
  double vzLo=-1e9, vzHi=1e9, icLo=-1, icHi=99999;
  double kcDenom=0, kcPivot=0, keOff=0, ebzSlope=0, ebzPivot=500;
  bool kcOn=false, ebzOn=false;
  int npLo=0, npHi=9999;
};
inline double omega2(double x,double y,double z){
  return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
inline double kcSlopeDeg(const Pars&p){
  return (p.kcDenom!=0 && std::fabs(p.kcDenom)>1e-9) ? 360.0/p.kcDenom : 0.0; }
inline double thCorr(const Pars&p,double th,double ke){
  return p.kcOn ? th - kcSlopeDeg(p)*(ke - p.kcPivot) : th; }
inline double ebAt(const Pars&p,double vz,bool hasVz){
  if(!p.ebzOn || !hasVz || !std::isfinite(vz)) return p.ebeam;
  return p.ebeam - (p.ebzSlope/1000.0)*(vz - p.ebzPivot); }
// returns Ex (MeV), or NaN below threshold
inline double exOf(const Pars&p,double Kp,double thlDeg,double Ke){
  const double M1=p.mBeam*U_AMU,M2=p.mTarg*U_AMU,M3=p.mEject*U_AMU,M4=p.mResid*U_AMU;
  double thl=thlDeg*M_PI/180.0;
  double Et1=Kp+M1, Et3=Ke+M3;
  double s=M1*M1+M2*M2+2*M2*Et1;
  double uu=M2*M2+M3*M3-2*M2*Et3;
  double a=(std::cos(thl)*omega2(s,M1*M1,M2*M2)*omega2(uu,M2*M2,M3*M3)
            -(s-M1*M1-M2*M2)*(M2*M2+M3*M3-uu))/(2*M2*M2)+s+uu-M2*M2;
  if(a<0) return std::nan("");
  return std::sqrt(a)-M4;
}
