#pragma once
#include "exlib.h"
// theta_CM of the detected ejectile, by boosting its lab 4-vector into the CM frame.
// For elastic recoil off a target at rest this reduces to theta_cm = 180 - 2*theta_lab, so the
// AT-TPC's backward lab angles (40-86 deg) map onto theta_cm ~ 8-100 deg -- the same domain the
// 15C(p,p') analysis bins over.
inline double thcmOf(const Pars& p, double Kp, double thlDeg, double Ke)
{
   const double M1 = p.mBeam * U_AMU, M2 = p.mTarg * U_AMU, M3 = p.mEject * U_AMU;
   double p1 = std::sqrt(Kp * Kp + 2 * Kp * M1), E1 = Kp + M1;
   TLorentzVector W(0, 0, p1, E1 + M2);
   double p3 = std::sqrt(Ke * Ke + 2 * Ke * M3), E3 = Ke + M3;
   double th = thlDeg * M_PI / 180.0;
   TLorentzVector d(p3 * std::sin(th), 0, p3 * std::cos(th), E3);
   d.Boost(-W.BoostVector());
   return d.Theta() * 180.0 / M_PI;
}
// The 15C(p,p') analysis bins the CM angle of the PROJECTILE; thcmOf gives the detected recoil,
// and for two-body kinematics the two are supplementary. Use this to match their domain.
inline double thcmProj(const Pars& p, double Kp, double thlDeg, double Ke)
{ return 180.0 - thcmOf(p, Kp, thlDeg, Ke); }
