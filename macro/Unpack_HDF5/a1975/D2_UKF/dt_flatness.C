/// \file dt_flatness.C
/// Ebeam-flatness metric for the a1975 16C(d,t)15C drift-velocity scan.
///
/// WHAT IT MEASURES.  For a trial beam energy the excitation energy is recomputed from the
/// CACHED (ke, theta_lab) of every triton -- the cache's own `ex` column is only valid at the
/// Ebeam it was built with, so it is deliberately NOT used here. Tracks near the isolated
/// 3.103 MeV level of 15C are selected, binned in theta, and the mean Ex of each bin is taken.
/// The metric is the spread (population rms) of those bin means:
///
///     rms = 0  <=>  the level sits at the same Ex at every angle  <=>  the kinematics close
///
/// The best Ebeam is the one that minimises it. This is a FLATNESS test, not an accuracy test:
/// the rms says the level does not drift with angle, it does NOT say the level landed at 3.103.
/// That is why `meanEx` is reported alongside -- at dv 1.50 the flattest Ebeam put the mean at
/// +4.005, i.e. the scan latched onto the wrong structure and the small rms meant nothing.
/// Always read rms and meanEx together.
///
/// WHY IT CAN RETURN NOTHING.  Every theta bin must hold at least `minPerBin` tracks. When the
/// Ex scale is stretched far enough (dv 1.60) the level is smeared out of the +-exWin window
/// and no Ebeam in the scan fills all the bins -- the macro then reports NONE rather than
/// inventing a minimum from two surviving bins. A "no result" is a real result here.
///
/// The kinematics below are copied VERBATIM from ex_dt_a1975.C (which took them from
/// C16_pp_ana.C) so that a cache and its scan cannot drift apart. `selfcheck` recomputes Ex at
/// the cache's own build energy and compares with the stored column; it must come back at the
/// 1e-3 MeV level, otherwise the copy has diverged or the cache was built at another Ebeam.
///
///   root -l -b -q 'dt_flatness.C("/mnt/f/a1975/caches/dt_kin_dv1136.root")'
///   root -l -b -q 'dt_flatness.C("...root",150,225,1,3.103,1.3,5,20,70,true,184)'

#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

static double omega2_f(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

// two-body kinematics (verbatim from ex_dt_a1975.C): returns {Ex, theta_cm[deg]}
static std::tuple<double, double> kine_2b_f(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                            double K_eject)
{
   double Et1 = K_proj + m1, Et2 = m2, Et3 = K_eject + m3, Et4 = Et1 + Et2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * omega2_f(s, m1 * m1, m2 * m2) * omega2_f(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg < 0)
      return {std::nan(""), std::nan("")};
   double m4_ex = std::sqrt(arg);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (omega2_f(s, m1 * m1, m2 * m2) * omega2_f(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

/// \param cache      per-dv kinematics ntuple built by cache_dt.sh
/// \param exRef      the level used as the ruler (3.103 MeV, isolated -- the 0/0.740 doublet is
///                   degenerate with any width being measured, so it is useless as a ruler)
/// \param exWin      half-window around exRef that defines "on the level"
/// \param nTh        number of theta bins the flatness is measured across
/// \param thLo,thHi  theta range those bins span
/// \param useThcm    bin in theta_cm (default) or theta_lab
/// \param cacheEbeam the energy the cache was BUILT at, used only by the self-check
void dt_flatness(TString cache, double ebLo = 150, double ebHi = 225, double ebStep = 1.0, double exRef = 3.103,
                 double exWin = 1.3, int nTh = 5, double thLo = 20, double thHi = 70, bool useThcm = true,
                 double cacheEbeam = 184.0, double chi2Max = 5, double icLo = 900, double icHi = 1300,
                 double keMin = 5, double vzLo = 50, double vzHi = 700, int minPerBin = 20, bool verbose = false)
{
   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;
   const double m_t = 3.01550072 * u, m_C15 = 15.0105993 * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", cache.Data());
      return;
   }
   auto *nt = dynamic_cast<TNtuple *>(f->Get("pk"));
   if (!nt) {
      printf("no ntuple 'pk' in %s\n", cache.Data());
      return;
   }
   float ke, theta, vertexz, vertexr, thcm, ex, chi2ndf, brho, dedx, sqrtdedx, ncl, ntrk, run, ic, kefit, thetafit,
      exfit;
   nt->SetBranchAddress("ke", &ke);
   nt->SetBranchAddress("theta", &theta);
   nt->SetBranchAddress("vertexz", &vertexz);
   nt->SetBranchAddress("ex", &ex);
   nt->SetBranchAddress("chi2ndf", &chi2ndf);
   nt->SetBranchAddress("ic", &ic);

   // load once: the Ebeam loop below re-evaluates the kinematics on these same tracks
   std::vector<double> vKe, vTh, vExStored;
   Long64_t nAll = nt->GetEntries();
   for (Long64_t i = 0; i < nAll; ++i) {
      nt->GetEntry(i);
      if (!(chi2ndf < chi2Max))
         continue;
      if (!(ic >= icLo && ic <= icHi))
         continue; // beam gate: off-cocktail events reconstruct ~1.7 MeV high
      if (!(ke > keMin))
         continue; // drops the low-KE pile-up at the kinematic turnover
      if (!(vertexz > vzLo && vertexz < vzHi))
         continue; // vertex fiducial
      vKe.push_back(ke);
      vTh.push_back(theta);
      vExStored.push_back(ex);
   }
   printf("\n%s\n  tracks: %lld in cache -> %zu after cuts "
          "(chi2ndf<%.3g, ic[%.0f,%.0f], ke>%.3g, vz[%.0f,%.0f])\n",
          gSystem->BaseName(cache.Data()), nAll, vKe.size(), chi2Max, icLo, icHi, keMin, vzLo, vzHi);
   if (vKe.empty())
      return;

   // --- self-check: reproduce the stored Ex at the cache's own build energy -------------------
   double d2 = 0;
   int nchk = 0;
   for (size_t i = 0; i < vKe.size(); ++i) {
      auto [exC, tcm] = kine_2b_f(m_C16, m_d, m_t, m_C15, cacheEbeam, vTh[i] * TMath::DegToRad(), vKe[i]);
      if (std::isnan(exC))
         continue;
      d2 += (exC - vExStored[i]) * (exC - vExStored[i]);
      ++nchk;
   }
   double selfrms = nchk ? std::sqrt(d2 / nchk) : -1;
   printf("  self-check vs stored ex at Ebeam=%.1f : rms %.4f MeV on %d tracks%s\n", cacheEbeam, selfrms, nchk,
          selfrms > 0.01 ? "   <-- WARNING: cache built at a different Ebeam?" : "");

   // --- Ebeam scan ---------------------------------------------------------------------------
   double bestEb = -1, bestRms = 1e9, bestMean = 0;
   int bestN = 0;
   int nValid = 0;
   printf("  scanning Ebeam %.0f..%.0f step %.2f | %d theta%s bins over [%.0f,%.0f], "
          "|Ex-%.3f|<%.2f, >=%d per bin\n",
          ebLo, ebHi, ebStep, nTh, useThcm ? "_cm" : "_lab", thLo, thHi, exRef, exWin, minPerBin);

   for (double eb = ebLo; eb <= ebHi + 1e-9; eb += ebStep) {
      std::vector<double> sum(nTh, 0.0), cnt(nTh, 0.0);
      double allSum = 0;
      int allN = 0;
      for (size_t i = 0; i < vKe.size(); ++i) {
         auto [exT, tcm] = kine_2b_f(m_C16, m_d, m_t, m_C15, eb, vTh[i] * TMath::DegToRad(), vKe[i]);
         if (std::isnan(exT))
            continue;
         if (std::fabs(exT - exRef) > exWin)
            continue;
         double th = useThcm ? tcm : vTh[i];
         int b = (int)((th - thLo) / (thHi - thLo) * nTh);
         if (b < 0 || b >= nTh)
            continue;
         sum[b] += exT;
         cnt[b] += 1;
         allSum += exT;
         ++allN;
      }
      bool full = true;
      for (int b = 0; b < nTh; ++b)
         if (cnt[b] < minPerBin)
            full = false;
      if (!full)
         continue; // not enough structure at this Ebeam to judge flatness
      ++nValid;
      double m = 0;
      for (int b = 0; b < nTh; ++b)
         m += sum[b] / cnt[b];
      m /= nTh;
      double v = 0;
      for (int b = 0; b < nTh; ++b) {
         double d = sum[b] / cnt[b] - m;
         v += d * d;
      }
      double rms = std::sqrt(v / nTh);
      if (verbose)
         printf("    Ebeam %6.1f  rms %.4f  meanEx %+7.3f  N %d\n", eb, rms, allSum / allN, allN);
      if (rms < bestRms) {
         bestRms = rms;
         bestEb = eb;
         bestMean = allSum / allN;
         bestN = allN;
      }
   }

   if (bestEb < 0) {
      printf("  ==> NONE: no Ebeam in the scan fills all %d theta bins with >=%d tracks.\n"
             "      The level structure is gone at this dv, not merely displaced.\n",
             nTh, minPerBin);
      return;
   }
   printf("  ==> best Ebeam %.1f | rms %.3f | meanEx %+.3f | N %d | %d/%d scan points valid%s\n", bestEb, bestRms,
          bestMean, bestN, nValid, (int)((ebHi - ebLo) / ebStep) + 1,
          std::fabs(bestMean - exRef) > 0.5 ? "\n      WARNING: meanEx is far from the reference level -- the scan is"
                                              " tracking the WRONG structure, the small rms is meaningless."
                                            : "");
}
