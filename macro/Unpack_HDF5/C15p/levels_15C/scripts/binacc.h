#pragma once
#include "excm.h"
// EXACT acceptance per (theta_CM bin, Ex).
// For two-body kinematics theta_CM and Ex determine theta_lab and KE completely, so the
// acceptance of a theta_CM bin at a given Ex is just the fraction of that bin whose (theta_lab,
// KE) lands in the detector window.  No free parameter, no solid-angle integral, no simulation.
//
// This also folds in the theta_corr the explorer applies, so it is the acceptance of the ANALYSIS,
// not of an idealised detector.
struct Win { double thLo, thHi, keLo, keHi; };
// lab angle and KE of the ejectile for a given CM angle of the projectile
inline void labOf(const Pars& p, double Kp, double Ex, double thcmProjDeg, double& thl, double& ke)
{
   const double M1 = p.mBeam * U_AMU, M2 = p.mTarg * U_AMU, M3 = p.mEject * U_AMU;
   double M4 = p.mResid * U_AMU + Ex;
   double p1 = std::sqrt(Kp * Kp + 2 * Kp * M1), E1 = Kp + M1;
   TLorentzVector W(0, 0, p1, E1 + M2);
   double s = W.M2(), rs = std::sqrt(s);
   if (rs <= M3 + M4) { thl = -1; ke = -1; return; }
   double pcm = omega2(s, M3 * M3, M4 * M4) / (2 * rs), E3 = std::sqrt(pcm * pcm + M3 * M3);
   // recoil CM angle = 180 - projectile CM angle
   double tr = (180.0 - thcmProjDeg) * M_PI / 180.0;
   TLorentzVector d(pcm * std::sin(tr), 0, pcm * std::cos(tr), E3);
   d.Boost(W.BoostVector());
   thl = d.Theta() * 180.0 / M_PI; ke = d.E() - M3;
}
inline double accBin(const Pars& p, double Kp, double Ex, double cmLo, double cmHi, const Win& w,
                     bool applyThCorr = true, int N = 120)
{
   int ok = 0;
   for (int i = 0; i < N; i++) {
      double tcm = cmLo + (cmHi - cmLo) * (i + 0.5) / N;
      double thl, ke; labOf(p, Kp, Ex, tcm, thl, ke);
      if (thl < 0) continue;
      double thc = applyThCorr ? thCorr(p, thl, ke) : thl;   // the analysis' own theta correction
      if (thc >= w.thLo && thc <= w.thHi && ke >= w.keLo && ke <= w.keHi) ok++;
   }
   return (double)ok / N;
}
