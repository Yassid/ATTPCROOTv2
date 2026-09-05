/// @file ex_genfit_3Hed.C
/// @brief Excitation-energy resolution of 46Ar(3He,d)47K from the GENFIT fit, by theta_lab band.
///
/// WHY THIS EXISTS. `compare_configs.C`, which produced the field x pad matrix quoted through
/// August, reads <tag>_reco.root and runs AtTools::AtSpyralPID -- the PRE-FIT first-arc circle
/// estimate. No genfit fit entered that table at all. Its headline conclusions, "2 mm pads yes,
/// 3.8 T no" and "the 3.8 T backward residual goes bimodal", are therefore statements about a
/// circle estimator, and the standing open item was to test them against a real fit. This macro
/// is that test: same kinematics, same binning, same median/IQR statistics, genfit input.
///
/// READS GetKinematicsXtr, NOT GetKinematics. With SetBackExtrapToAxis on -- the default in
/// fitGenfitter_Ar46.C -- the vertex-gap-corrected momentum lands in the Xtr slot and the raw
/// genfit momentum stays in the other one. Reading the wrong slot silently loses the correction.
/// Pass useXtr = kFALSE to read the raw slot instead; the macro prints which it used.
///
/// MEDIAN AND IQR, NOT MEAN AND FWHM. These distributions are a narrow core on broad tails. A
/// half-maximum walk out from the peak returned 0.765 against 1.800 MeV for two histograms that
/// differed only by a CONSTANT shift of every entry, and the earlier "matEffects halves the width"
/// claim came from exactly that. Do not reintroduce it.
///
/// DRIFT-Z IS MIRRORED in this simulation (r = -1.000 against truth, z_true + z_reco = 101.2 cm),
/// so the beam energy at the vertex needs zUse = driftLength - z_reco. The sign is not assumed:
/// it is measured per configuration from the correlation and printed.
///
/// The defaults are the CURRENT three-field ground-state comparison (2.0 / 2.85 / 3.8 T on real
/// per-field Magboltz transport). Note simDirs: the 2.85 and 3.8 T arms re-reconstruct from sims
/// that live in the OLD placeholder directories, because generation depends on the B field and
/// never on the .par, so those sims stayed valid when the transport was corrected.
///
///   root -b -q 'ex_genfit_3Hed.C()'                                  // the adopted comparison
#include "ex_core_3Hed.h"   // the shared inversion -- see that file before changing anything here

void ex_genfit_3Hed(TString tags = "gs_s3001,gs_s3002",
                    TString dirs = "/mnt/f/ar46_3hed_mb_B20,/mnt/f/ar46_3hed_mb_B285,/mnt/f/ar46_3hed_mb_B38",
                    TString names = "2.0 T,2.85 T,3.8 T",
                    Double_t dThetaMax = 10.0, Double_t driftLength = 100.0, Bool_t useXtr = kTRUE,
                    Double_t chi2Max = -1.0,
                    // The _mb_ arms at 2.85 and 3.8 T re-reconstruct from sims that live in the
                    // OLD directories (generation depends on the field, not the .par, so it is
                    // reused). Their _sim.root is therefore NOT next to their fit. Empty = same
                    // directory as the fit.
                    TString simDirs = "/mnt/f/ar46_3hed_mb_B20,/mnt/f/ar46_3hed_OLD_2.85T_placeholder,"
                                      "/mnt/f/ar46_3hed_OLD_3.8T_placeholder")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");

   TObjArray *aT = tags.Tokenize(","), *aD = dirs.Tokenize(","), *aN = names.Tokenize(",");
   TObjArray *aS = simDirs.Length() ? simDirs.Tokenize(",") : nullptr;
   const int NB = 8;
   double edge[NB + 1] = {50, 60, 70, 80, 90, 100, 110, 125, 145};

   printf("\nExcitation energy from GENFIT, %s slot%s\n", useXtr ? "GetKinematicsXtr (back-extrapolated)" : "GetKinematics (raw)",
          chi2Max > 0 ? Form(", chi2/ndf < %.2f", chi2Max) : ", no chi2 cut");

   for (int ic = 0; ic < aD->GetEntries(); ++ic) {
      TString dir = ((TObjString *)aD->At(ic))->GetString();
      TString nm = (ic < aN->GetEntries()) ? ((TObjString *)aN->At(ic))->GetString() : dir;

      TString sdir = (aS && ic < aS->GetEntries()) ? ((TObjString *)aS->At(ic))->GetString() : dir;
      Ar46::Sample S = Ar46::Collect(dir, sdir, aT, useXtr, dThetaMax, chi2Max, driftLength, nm);

      std::vector<double> band[NB], all = S.ex;
      for (size_t j = 0; j < S.ex.size(); ++j)
         for (int b = 0; b < NB; ++b)
            if (S.thetaTrue[j] >= edge[b] && S.thetaTrue[j] < edge[b + 1]) band[b].push_back(S.ex[j]);
      const long nTruth = S.nTruth, nFit = S.nFit, nCut = S.nCut;
      const double rV = S.rV;
      const bool useMirror = S.mirror;

      printf("\n===== %s  (%s) =====\n", nm.Data(), dir.Data());
      printf("  truth deuterons %ld | matched fits %ld (%.1f%%)%s | drift-z mirror %s\n", nTruth, nFit,
             100.0 * nFit / std::max(1L, nTruth),
             chi2Max > 0 ? Form(" | chi2-cut %ld", nCut) : "", useMirror ? "APPLIED" : "not applied");
      printf("  vertex-z correlation against truth r = %+.3f  (negative => the sim mirrors drift z)\n", rV);
      printf("  %-12s %7s %10s %10s\n", "theta_lab", "n", "med Ex", "IQR");
      for (int b = 0; b < NB; ++b) {
         if (band[b].size() < 8) continue;
         printf("  %4.0f-%-7.0f %7zu %10.3f %10.3f\n", edge[b], edge[b + 1], band[b].size(), Ar46::MedOf(band[b]),
                Ar46::IqrOf(band[b]));
      }
      printf("  %-12s %7zu %10.3f %10.3f\n", "ALL", all.size(), Ar46::MedOf(all), Ar46::IqrOf(all));
   }
   printf("\nNOTE: these tags mix the three states (0, 0.36, 2.02 MeV), so med Ex is a mixture\n"
          "centroid and only the IQR is a resolution. Run one state's tags alone for a centroid.\n");
}
