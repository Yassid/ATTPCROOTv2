#pragma once
#include "excm.h"
#include "binacc.h"
// =============================================================================================
//  Is the LOW-ENERGY BACKGROUND fusion evaporation?
//
//  ⚠ What has already been tested and REJECTED (spyral-15Cd-calibration): evaporation as the
//  explanation for the FLAT 2-20 MeV CONTINUUM -- that spectrum is flat where a Maxwellian is
//  not, chi2/ndf = 1889.  The component tested here is a DIFFERENT one: the background that turns
//  on near 2.5 MeV, decays with kappa ~ 0.24-0.30 /MeV, sits in the FORWARD theta_CM bins and is
//  made of real deuterons at KE ~ 1.4 MeV.  That has never been tested against evaporation.
//
//  ⚠ Two tests that look decisive and are NOT (same source, do not repeat them):
//    (1) "fraction below the evaporation kinematic limit" -- the limit lies exactly ON the elastic
//        locus because 17N* -> 15C(gs) + d is two-body, so it discriminates nothing;
//    (2) a compound-frame isotropy test on a BOX in (KE, theta) -- the box forces the answer.
//  Use the ENERGY SPECTRUM SHAPE over the full sample, which is what this does.
//
//  MODEL.  15C + d fuses to 17N* which evaporates a deuteron isotropically in its own frame with
//  a surface-Maxwellian spectrum above the Coulomb barrier,
//        dN/deps  ~  (eps - B) * exp(-(eps - B)/T),     eps > B
//  then the deuteron is boosted to the lab, passed through the analysis' own cuts, and its Ex is
//  reconstructed with the SAME two-body deuteron hypothesis the data uses -- so the simulated and
//  measured quantities are the same thing.
namespace EV {
   const double U = 931.49401;
   const double M15C = 15.0105993*U, MD = 2.0141018*U, M17N = 17.008449*U;
   // compound excitation energy and temperature
   inline double Estar(double Kb){
      double pb=std::sqrt(Kb*Kb+2*Kb*M15C);
      TLorentzVector W(0,0,pb,Kb+M15C+MD);
      double Ecm=W.M()-M15C-MD;                        // CM kinetic energy
      double Q = M15C+MD-M17N;                         // fusion Q value
      return Ecm+Q; }
   inline double temperature(double Estar,double A=17.0){
      double a=A/8.0; return std::sqrt(std::max(0.0,Estar)/a); }
   inline double coulombBarrier(){
      double R=1.4*std::pow(17.0,1.0/3.0)+2.0;         // fm, CN radius + deuteron
      return 1.0*7.0*1.44/R; }                          // Z_d * Z_(N) * e^2/R, MeV
}
