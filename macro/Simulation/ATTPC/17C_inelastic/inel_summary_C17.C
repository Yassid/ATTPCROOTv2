/// @file inel_summary_C17.C
/// @brief The 17C(p,p') / 17C(d,d') result: excitation-energy resolution, the method floor beside
///        it, the vertex beam-energy correction, and whether the 1/2+ / 5/2+ doublet separates --
///        across the campaign's four cells (two probes x two fields).
///
///   root -b -q 'inel_summary_C17.C'
///   root -b -q 'inel_summary_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel")'
///
/// WHAT IT REPORTS, and why each is here rather than a single number:
///
///  1. sigma(Ex) as IQR/1.349 -- a robust width, not an RMS and not a fit. Fit-failure tails in
///     this pipeline are heavy enough to fake large effects and a median-based width ignores them.
///  2. THE METHOD FLOOR beside every measured width: the same inversion applied to the TRUTH
///     (keTrue, thTrue) at the same constant beam energy. floor ~ measured means the detector is
///     invisible and no detector change helps; floor << measured means there is something to
///     chase. A resolution quoted without it cannot tell those apart.
///  3. THE VERTEX BEAM-ENERGY CORRECTION. The 17C beam loses ~14.5 MeV crossing the chamber, so one
///     constant Ebeam is wrong by up to +-7 MeV in a way perfectly correlated with a quantity the
///     detector already measures. E_beam(z) is MEASURED here from truth, not assumed.
///  4. LEVEL SEPARATION as Delta(Ex)/(sigma_1 + sigma_2), for g.s.->217 (217 keV) and the pair the
///     proposal lives on, 217->332 (115 keV). Above ~2 a pair is resolvable in a fit; below ~1 it
///     is one bump. This macro does NOT decide whether the proposal works -- at separation < 1 the
///     question becomes whether a FIXED-position fit can still recover the two amplitudes, which
///     is decompose_C17.C.
///
/// THE REFERENCE IS THE GENERATED Ex, NOT THE TREE'S exTrue. ex_res_C14_hf.C:162 builds exTrue by
/// putting the truth kinematics through this same inversion at this same constant beam energy, so
/// solving anything against it is self-fulfilling -- it would return the constant exactly and a
/// floor measured that way is identically zero. Both happened while building the (d,p) arm.
/// exTrue is then free: it already IS the method-floor residual.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
double s17_omega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

/// Two-body inversion, identical in form to acceptance_C14.C:acc_kine -- given the measured
/// ejectile (theta_lab, KE) and a beam energy, return the residual excitation energy.
double s17_ex(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * s17_omega2(s, m1 * m1, m2 * m2) * s17_omega2(u, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                         (2 * m2 * m2) +
                      s + u - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

/// robust width: interquartile range over 1.349, equal to sigma for a gaussian but immune to the
/// fit-failure tail. -1 on too little data.
double s17_iqrSigma(std::vector<double> v)
{
   if (v.size() < 20)
      return -1;
   std::sort(v.begin(), v.end());
   return (v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())]) / 1.349;
}
double s17_median(std::vector<double> v)
{
   if (v.empty())
      return 0;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}
} // namespace

void inel_summary_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", Double_t EbeamConst = 135.0,
                      Double_t chi2Cut = 5.0)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double mBeam = 17.0225787 * u; // 17C
   const double mRes = mBeam;           // scattering: the residual IS the beam

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};
   const char *stName[nL] = {"3/2+ g.s.", "1/2+ 217", "5/2+ 332"};

   struct Cfg {
      const char *tag;
      const char *label;
      double mLight; // target = ejectile, amu
   };
   const int nC = 4;
   const Cfg cfg[nC] = {{"pp_b285", "17C(p,p')  2.85 T", 1.007825},
                        {"pp_b400", "17C(p,p')  4.00 T", 1.007825},
                        {"dd_b285", "17C(d,d')  2.85 T", 2.0141018},
                        {"dd_b400", "17C(d,d')  4.00 T", 2.0141018}};

   // The light recoil of an inelastic channel cannot exceed 90 deg in the lab, so these slices are
   // the physical range and not an inherited restriction (unlike in the (d,p) arm).
   const int nSl = 7;
   const double slLo[nSl] = {20, 30, 40, 50, 60, 70, 80};
   const double slHi[nSl] = {30, 40, 50, 60, 70, 80, 90};

   printf("\n\033[1;33m##########################################################################\033[0m\n");
   printf("\033[1;33m 17C(p,p') and 17C(d,d') -- 300 torr H2 / D2, real AT-TPC pad plane\033[0m\n");
   printf("\033[1;33m constant-Ebeam analysis at %.2f MeV, chi2/ndf < %.1f\033[0m\n", EbeamConst, chi2Cut);
   printf("\033[1;33m##########################################################################\033[0m\n");
   printf("  the pair that decides the proposal is 217 -> 332, a %.0f keV gap\n", (ExGen[2] - ExGen[1]) * 1000);

   // collected for the cross-cell table at the end
   double sigC[nC][nL], sigV[nC][nL], sigF[nC][nL];
   double sepC[nC], sepV[nC];
   for (int c = 0; c < nC; ++c) {
      for (int l = 0; l < nL; ++l)
         sigC[c][l] = sigV[c][l] = sigF[c][l] = -1;
      sepC[c] = sepV[c] = -1;
   }

   for (int c = 0; c < nC; ++c) {
      const double m2 = cfg[c].mLight * u; // target
      const double m3 = m2;                // ejectile
      printf("\n\033[1;36m================ %s ================\033[0m\n", cfg[c].label);

      std::vector<double> resC[nL], resV[nL], floorC[nL];
      std::vector<double> resCsl[nL][nSl], resVsl[nL][nSl], floorCsl[nL][nSl];
      // tracking residuals, for the closure test against the kinematics gate
      std::vector<double> dKE[nL], dTh[nL];
      std::vector<double> dKEsl[nL][nSl];
      long nTot[nL] = {0}, nTail[nL] = {0};
      double ebz_a = 0, ebz_b = 0;
      bool haveEbz = false;

      for (int l = 0; l < nL; ++l) {
         // the tag written by ex_res_C14_hf.C is the full job name: <chan>_<state>_<btag>_s<seed>
         TString chan(cfg[c].tag);
         chan.Remove(2); // "pp" or "dd"
         TString btag(cfg[c].tag);
         btag.Remove(0, 3); // "b285" or "b400"
         // Quote ONLY the directory. The products live under ".../Seagate Hub/...", so the path
         // needs quoting for the space -- but quoting the whole pattern stops the shell expanding
         // the "s*" wildcard and every file silently reads as MISSING.
         TString dirq = TString("\"") + root + "/" + cfg[c].tag + "\"";
         TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
         TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
         found = found.Strip(TString::kBoth);
         if (found.IsNull()) {
            printf("\033[1;31m  %-10s MISSING (%s)\033[0m\n", stTag[l], pat.Data());
            continue;
         }
         TFile *f = TFile::Open(found);
         if (!f || f->IsZombie()) {
            printf("\033[1;31m  cannot open %s\033[0m\n", found.Data());
            continue;
         }
         TTree *t = (TTree *)f->Get("res");
         if (!t) {
            printf("\033[1;31m  no res tree in %s\033[0m\n", found.Data());
            continue;
         }
         double exR, exT, thT, thR, keT, keR, cmT, c2n, zT, zR;
         t->SetBranchAddress("exReco", &exR);
         t->SetBranchAddress("exTrue", &exT);
         t->SetBranchAddress("thTrue", &thT);
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keTrue", &keT);
         t->SetBranchAddress("keReco", &keR);
         t->SetBranchAddress("cmTrue", &cmT);
         t->SetBranchAddress("chi2ndf", &c2n);
         t->SetBranchAddress("zTrue", &zT);
         t->SetBranchAddress("zReco", &zR);

         // ---- MEASURE E_beam(z) from truth on the elastic sample. For each event, solve by
         // ---- bisection for the beam energy that makes the GENERATED Ex close exactly, then fit
         // ---- against the true vertex z. Nothing is taken from a stopping-power table.
         // Any level will do: the bisection solves for the beam energy that reproduces THAT level's
         // generated Ex, so it does not have to be the elastic. Keyed on the first level actually
         // present, so a partially-finished campaign still yields the correction.
         if (!haveEbz) {
            TGraph g;
            const double tgt = ExGen[l];
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               double lo = EbeamConst - 40, hi = EbeamConst + 40;
               if ((s17_ex(mBeam, m2, m3, mRes, lo, thT, keT) - tgt) *
                      (s17_ex(mBeam, m2, m3, mRes, hi, thT, keT) - tgt) >
                   0)
                  continue;
               for (int it = 0; it < 60; ++it) {
                  const double mid = 0.5 * (lo + hi);
                  if ((s17_ex(mBeam, m2, m3, mRes, lo, thT, keT) - tgt) *
                         (s17_ex(mBeam, m2, m3, mRes, mid, thT, keT) - tgt) <=
                      0)
                     hi = mid;
                  else
                     lo = mid;
               }
               g.SetPoint(g.GetN(), zT, 0.5 * (lo + hi));
            }
            if (g.GetN() > 100) {
               TF1 lin("lin", "pol1");
               g.Fit(&lin, "QN");
               ebz_a = lin.GetParameter(0);
               ebz_b = lin.GetParameter(1);
               haveEbz = true;
               printf("  E_beam(z) measured from truth: %.4f %+.6f * z[mm]  (n = %d)\n", ebz_a, ebz_b, g.GetN());
               printf("    -> %.2f MeV at z = 0, %.2f at z = 1000; loss over the metre %.2f MeV "
                      "(the constant used is %.2f)\n",
                      ebz_a, ebz_a + 1000 * ebz_b, -1000 * ebz_b, EbeamConst);
            }
         }

         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (c2n >= chi2Cut)
               continue;
            ++nTot[l];
            // residual against the GENERATED excitation. exR is built with the GROUND-STATE
            // residual mass -- an experiment does not know the level in advance -- so it should
            // come out at ExGen and this difference is the error.
            const double rC = exR - ExGen[l];
            const double fC = exT; // already the method-floor residual; see the header note
            double rV = rC;
            if (haveEbz)
               rV = s17_ex(mBeam, m2, m3, mRes, ebz_a + ebz_b * zR, thR, keR) - ExGen[l];
            resC[l].push_back(rC);
            resV[l].push_back(rV);
            floorC[l].push_back(fC);
            dKE[l].push_back(keR - keT);
            dTh[l].push_back(thR - thT);
            if (std::fabs(rC) > 1.0)
               ++nTail[l];
            for (int s = 0; s < nSl; ++s)
               if (thR >= slLo[s] && thR < slHi[s]) {
                  resCsl[l][s].push_back(rC);
                  resVsl[l][s].push_back(rV);
                  floorCsl[l][s].push_back(fC);
                  dKEsl[l][s].push_back(keR - keT);
               }
         }
         f->Close();
      }

      // ---- per-level widths -------------------------------------------------------------------
      printf("\n  %-12s %8s %10s %10s %10s %9s %9s\n", "level", "n", "sigma_C", "sigma_V", "floor", "median_C",
             "tail>1MeV");
      for (int l = 0; l < nL; ++l) {
         if (!nTot[l])
            continue;
         sigC[c][l] = s17_iqrSigma(resC[l]);
         sigV[c][l] = s17_iqrSigma(resV[l]);
         sigF[c][l] = s17_iqrSigma(floorC[l]);
         printf("  %-12s %8ld %10.3f %10.3f %10.3f %9.3f %8.1f %%\n", stName[l], nTot[l], sigC[c][l], sigV[c][l],
                sigF[c][l], s17_median(resC[l]), 100.0 * nTail[l] / nTot[l]);
      }

      // ---- closure against the kinematics gate --------------------------------------------------
      // inel_kinematics_C17.C propagated the 14C(p,p') tracking resolution through this channel's
      // leverage before any event existed, and predicted a sigma(Ex) FLAT IN ANGLE. The leverage
      // half of that is right -- |dEx/dKE| really is flat at 0.533 (p,p') / 0.563 (d,d') -- but the
      // flatness conclusion was WRONG, because it held sigma(KE) at the single elastic-average
      // number the 14C matrix quotes. sigma(KE) is strongly angle-dependent: at large theta_lab the
      // recoil is slow, its helix is tight, and it is measured far better. That is why the table
      // below is printed per slice, and why the good resolution lives at large theta_lab (= small
      // theta_cm) rather than everywhere.
      {
         const double lev = (cfg[c].mLight < 1.5) ? 0.533 : 0.563;
         const double sKE = s17_iqrSigma(dKE[0]), sTh = s17_iqrSigma(dTh[0]);
         printf("\n  tracking (elastic, all angles): sigma(KE_reco - KE_true) = %.3f MeV, "
                "sigma(theta) = %.3f deg\n",
                sKE, sTh);
         printf("  %-14s", "slice");
         for (int s = 0; s < nSl; ++s)
            printf(" %10.0f-%.0f", slLo[s], slHi[s]);
         printf("\n");
         for (int l = 0; l < nL; ++l) {
            if (!nTot[l])
               continue;
            printf("  %-9s %-4s", l == 0 ? "sigma(KE)" : "", stTag[l]);
            for (int s = 0; s < nSl; ++s) {
               const double v = s17_iqrSigma(dKEsl[l][s]);
               if (v < 0)
                  printf(" %13s", "-");
               else
                  printf(" %13.3f", v);
            }
            printf("\n");
         }
         printf("  %-14s", "-> x |dEx/dKE|");
         for (int s = 0; s < nSl; ++s) {
            const double v = s17_iqrSigma(dKEsl[0][s]);
            if (v < 0)
               printf(" %13s", "-");
            else
               printf(" %13.3f", lev * v);
         }
         printf("   <- predicted sigma(Ex), to compare with the F/V rows below\n");
      }

      // ---- per-slice, on the two proposal states ----------------------------------------------
      printf("\n  sigma(Ex) per theta_lab slice   [C = constant Ebeam, V = vertex-corrected, F = floor]\n");
      printf("  %-12s", "slice");
      for (int s = 0; s < nSl; ++s)
         printf(" %10.0f-%.0f", slLo[s], slHi[s]);
      printf("\n");
      for (int l = 0; l < nL; ++l) {
         if (!nTot[l])
            continue;
         for (int k = 0; k < 3; ++k) {
            printf("  %-9s %-2s", k == 0 ? stName[l] : "", k == 0 ? "C" : (k == 1 ? "V" : "F"));
            for (int s = 0; s < nSl; ++s) {
               const double v = k == 0 ? s17_iqrSigma(resCsl[l][s])
                                       : (k == 1 ? s17_iqrSigma(resVsl[l][s]) : s17_iqrSigma(floorCsl[l][s]));
               if (v < 0)
                  printf(" %13s", "-");
               else
                  printf(" %13.3f", v);
            }
            printf("\n");
         }
      }

      // ---- separation --------------------------------------------------------------------------
      printf("\n  level separation, Delta / (sigma_a + sigma_b)\n");
      printf("  %-22s %10s %14s %14s\n", "pair", "Delta[MeV]", "const Ebeam", "vertex-corr");
      for (int l = 0; l + 1 < nL; ++l) {
         if (sigC[c][l] < 0 || sigC[c][l + 1] < 0)
            continue;
         const double d = ExGen[l + 1] - ExGen[l];
         const double sc = d / (sigC[c][l] + sigC[c][l + 1]);
         const double sv = d / (sigV[c][l] + sigV[c][l + 1]);
         printf("  %-22s %10.3f %14.2f %14.2f\n", Form("%s -> %s", stTag[l], stTag[l + 1]), d, sc, sv);
         if (l == 1) {
            sepC[c] = sc;
            sepV[c] = sv;
         }
      }
   }

   // ---- the cross-cell table, which is the thing to quote --------------------------------------
   printf("\n\033[1;33m########## the campaign in one table ##########\033[0m\n");
   printf("  sigma(Ex) [MeV] and the separation of the 115 keV proposal doublet\n\n");
   printf("  %-20s %12s %12s %12s %14s %14s\n", "configuration", "sigma(gs)", "sigma(217)", "sigma(332)",
          "sep 217-332 C", "sep 217-332 V");
   for (int c = 0; c < nC; ++c) {
      printf("  %-20s", cfg[c].label);
      for (int l = 0; l < nL; ++l) {
         if (sigV[c][l] < 0)
            printf(" %12s", "-");
         else
            printf(" %12.3f", sigV[c][l]);
      }
      if (sepC[c] < 0)
         printf(" %14s %14s\n", "-", "-");
      else
         printf(" %14.2f %14.2f\n", sepC[c], sepV[c]);
   }
   printf("\n  (sigma columns are VERTEX-CORRECTED; the constant-Ebeam values are in the per-cell\n");
   printf("   tables above. Separation above ~2 is resolvable in a fit, below ~1 is one bump --\n");
   printf("   at which point decompose_C17.C is what decides whether the proposal still works.)\n");
   printf("\n  inel summary done\n\n");
}
