/// @file reco_dlc_optimal.C
/// @brief OPTIMAL PUMA reconstruction with the DLC resistive anode ON.
///
///  Winning chain (see PUMA_report.tex §7 / PUMA_dlc_findings.pdf):
///     AtPSAMultiFit (impulse-time z)      -- accurate per-pad charge (fitted integral)
///     -> HDBSCAN track finder             -- clusters the dense ~4x DLC band cleanly
///     -> per-ring charge-weighted centroid-- recovers the sub-pad position
///     -> UKF + GENFIT Kalman fits         -- GENFIT (0.3 mm meas.) exploits the sub-pad info
///     + Cu-trap/cryostat material corr.   -- removes the upstream energy-loss bias
///
///  Expected @375 MeV/c (pi+pi-):  GENFIT bias ~0%, sigma ~7% (2x better than no-DLC);
///                                 UKF sigma ~15%, bias ~-6% (calibratable).
///
///  This is a thin, documented wrapper over run_digi_ukf_genfit_test8.C that hard-sets
///  the winning 19-argument configuration so it is reproducible in one call.
///
///  Run: root -b -q 'reco_dlc_optimal.C(500,"data/attpcsim.root","data/dlc_optimal.root")'
void reco_dlc_optimal(int nEvents = 500, TString simFile = "data/attpcsim.root",
                      TString outFile = "data/dlc_optimal.root", TString species = "pi",
                      double rDLC = 1.35e6 /* 1.35 MOhm/sq, PUMA DLC sheet resistance */)
{
   gSystem->Setenv("PUMA_IN", simFile.Data());
   gROOT->LoadMacro("run_digi_ukf_genfit_test8.C");
   //  n , tCluster, gfFlip, PSA        , species, clTgt, nRings, nPads, prfSig, ring, primary, merge, minHits, backExtrapMat, matCorr, rDLC, cDLC   , thr, praFinder
   TString cmd = Form("run_digi_ukf_genfit_test8(%d, 8.0, false, \"mfimpulse\", \"%s\", 8, 16, 256, "
                      "0.0, true, true, false, 20, false, 18.0, %g, 6.0e-13, 0.0, \"hdbscan\")",
                      nEvents, species.Data(), rDLC);
   gROOT->ProcessLine(cmd);
   gSystem->Rename("data/output_digi_both8.root", outFile.Data());
   printf("\nDLC_OPTIMAL_DONE -> %s\n", outFile.Data());
   printf("Analyse with:  root -b -q 'compare_ukf_genfit_test8.C(\"%s\",<E_GeV>,\"./%s\",\"./%s\")'\n",
          species.Data(), outFile.Data(), simFile.Data());
}
