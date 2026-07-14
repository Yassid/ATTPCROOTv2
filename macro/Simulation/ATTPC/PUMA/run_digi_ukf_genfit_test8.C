/// @file run_digi_ukf_genfit_test8.C
/// @brief Digitize the PUMA branch-8 test sample (back-to-back pi+/pi-) and fit
///        it with BOTH Kalman filters in one pass for an apples-to-apples
///        comparison on identical clusters:
///          AtTrackingEventUKF     — EventFit::AtFitterUKF (unscented KF + CATIMA)
///          AtTrackingEventGenfit  — EventFit::AtGenfitter (genfit KalmanFitterRefTrack,
///                                   native Bethe-Bloch, per-track charge from PRA)
///        Both read the same AtPatternEvent, so their fitted tracks are directly
///        comparable against the known branch-8 truth (|p| = 0.3749 GeV/c, theta = 90 deg).
///
/// Chain: clusterize -> pulse -> PSA -> PRA -> UKF(pi+/pi-) -> genfit(pi+/pi-).
/// Reads ./data/attpcsim.root, writes ./data/output_digi_both8.root.
///
/// PSA default is "mfimpulse": AtPSAMultiFit fitting the same GET response the digi
/// used, with z from the fitted impulse time t0 — this removes the +8.6 mm shaping-
/// delay z bias AT THE SOURCE, so no SetVertexZBias correction is applied. Pass
/// psaType="max" to use AtPSAMax (peak) instead, which re-enables the +8.6 mm bias
/// correction for a back-to-back comparison.
///
/// Run: root -b -q 'run_digi_ukf_genfit_test8.C(1000)'                     // mfimpulse PSA
///      root -b -q 'run_digi_ukf_genfit_test8.C(1000, 8.0, false, "max")'  // AtPSAMax PSA
///      root -b -q 'run_digi_ukf_genfit_test8.C(1000, 8.0, true)'          // flip genfit charge-sign map

/// species: "pi" (branch 8) or "K"/"kaon" (branch 10). Sets the fit mass/PDG for
/// both fitters. Both charge hypotheses are offered to the UKF; the PRA charge
/// prior gates the sign, so a same-charge final state (branch 10 K+K+) is fine.
void run_digi_ukf_genfit_test8(int nEvents = 1000, float tCluster = 8.0, bool genfitChargeFlip = false,
                               TString psaType = "mfimpulse", TString species = "pi", int clusterTarget = 8,
                               int nRings = 16, int nPads = 256, double prfSigma = 0.0, bool ringClustering = false,
                               bool primaryOnly = false, bool useMerge = false, int minHits = 0,
                               bool backExtrapMat = false, double matCorrRefDp = 0.0,
                               double rDLC = 0.0, double cDLC = 6.0e-13, double psaThreshold = 0.0,
                               TString praFinder = "triplclust")
{
   const bool isK = (species == "K" || species == "kaon");
   // "both" registers K+/K-/pi+/pi- UKF hypotheses for a chi2-PID baseline vs the
   // dE/dx PID (pid_dedx.C). genfit + p0-mass fall back to pion in that mode.
   const bool isBoth = (species == "both");
   const double mMeV = isK ? 493.677 : 139.57039; // kaon / pion mass (pion for "both")
   const int pdgPos = isK ? 321 : 211;            // K+ / pi+
   const std::string nm = isK ? "K" : "pi";

   TString inOutDir = "./data/";
   // Output dir override (PUMA_OUT env): lets the sweep write ~GB digi files straight
   // to a big drive (e.g. F:) instead of growing the WSL vhdx on C:. Input still ./data.
   TString outDir = gSystem->Getenv("PUMA_OUT");
   if (outDir.IsNull()) outDir = inOutDir;
   TString outputFile = outDir + "output_digi_both8.root";
   TString paramFile = "ATTPC.PUMA_sim.par";

   TString dir = getenv("VMCWORKDIR");
   // Input file override (PUMA_IN env): point the digi at a specific sim file
   // (e.g. attpcsim_K.root) without editing the macro. Default ./data/attpcsim.root.
   TString mcFile = gSystem->Getenv("PUMA_IN");
   if (mcFile.IsNull()) mcFile = inOutDir + "attpcsim.root";
   TString digiParFile = dir + "/parameters/" + paramFile;

   // genfit's TGeoMaterialInterface needs gGeoManager (FAIRGeom) in memory.
   TString geoManFile = dir + "/geometry/ATTPC_PUMA_geomanager.root";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);
   fRun->SetGeomFile(geoManFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   if (gROOT->FindObject("FAIRGeom") == nullptr) { // genfit needs gGeoManager
      TFile *gf = TFile::Open(geoManFile);
      if (gf)
         gf->Get("FAIRGeom");
   }

   // PUMA annular pad plane: 16 equal-area rings x 256 pads (4096), R=62.9-121.1 mm.
   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, nRings, nPads);
   mapping->GeneratePadPlane();

   // primaryOnly: drop trackID<0 secondary/delta-ray deposits at clusterization so only
   // the primary tracks are digitized (clean pi arcs, no delta-ray junk tracks in PRA).
   auto clusterizeAlgo = std::make_shared<AtClusterize>();
   clusterizeAlgo->SetPrimaryOnly(primaryOnly);
   AtClusterizeTask *clusterizer = new AtClusterizeTask(clusterizeAlgo);
   clusterizer->SetPersistence(kFALSE);

   auto atPulse = std::make_shared<AtPulse>(mapping);
   if (rDLC > 0) {
      // DLC-physical path: derive the PRF sigma from the measured DLC sheet
      // resistance rDLC [Ohm/sq] (PUMA 1.2-1.5e6) and areal capacitance cDLC
      // [F/mm^2]; t=0.5 us ~ shaping time. Overrides prfSigma when set.
      atPulse->SetChargeDispersionFromDLC(rDLC, cDLC, 0.5);
      std::cout << "[digi] DLC dispersion from R=" << rDLC / 1e6 << " MOhm/sq" << std::endl;
   } else if (prfSigma > 0) {
      atPulse->SetChargeDispersion(prfSigma); // resistive-pad PRF (mm) — enables sub-pad centroiding
   }
   AtPulseTask *pulse = new AtPulseTask(atPulse);
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   // PSA: "max" = AtPSAMax (peak; carries the +8.6 mm shaping-delay z bias, corrected
   // downstream via SetVertexZBias). "mfimpulse" = AtPSAMultiFit with z from the fitted
   // IMPULSE time t0 -> removes the shaping-delay z offset AT THE SOURCE. NOTE: the fitted
   // t0 is poorly constrained for these fast pulses -> hit z-RMS ~15 mm (vs 1.3 mm for
   // AtPSAMax) -> PRA 3D clustering merges/fragments the back-to-back tracks -> momentum
   // resolution DOUBLES (20%->46%). "mfpeak" = same AtPSAMultiFit multi-peak deconvolution
   // but z from the fitted-curve PEAK (robust, like AtPSAMax) -> keeps multi-peak splitting
   // for overlapping pulses WITHOUT the noisy-z penalty.
   const bool useMFimpulse = (psaType == "mfimpulse");
   const bool useMFpeak = (psaType == "mfpeak");
   AtPSAtask *psaTask = nullptr;
   if (useMFimpulse || useMFpeak) {
      auto mf = std::make_unique<AtPSAMultiFit>();
      mf->SetThreshold(psaThreshold);
      mf->SetPeakingTime(0.5);              // 500 ns, matches PUMA AtPulse shaping (AtNominalResponse)
      mf->SetUseImpulseTime(useMFimpulse);  // impulse t0 (noisy) vs fitted peak (robust)
      psaTask = new AtPSAtask(std::move(mf));
   } else {
      auto psa = std::make_unique<AtPSAMax>();
      psa->SetThreshold(psaThreshold);
      psaTask = new AtPSAtask(std::move(psa));
   }
   psaTask->SetPersistence(kTRUE);

   AtPRAtask *praTask = nullptr;
   if (praFinder == "hdbscan") {
      // density clustering — separates dense DLC bands / same-charge pairs better
      auto finder = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
      finder->SetChargeFromCenter(true);
      praTask = new AtPRAtask(std::move(finder));
      praTask->SetMinNumHits(minHits > 0 ? minHits : 6);
   } else {
      praTask = new AtPRAtask(); // TriplClust (AtTrackFinderTC)
      praTask->SetTcluster(tCluster);
      praTask->SetTCUseSelectAndMerge(false); // blunt end-proximity merge over-fuses 2 tracks at 375
      praTask->SetTCUseCircleMerge(useMerge); // annular-safe: merge only SAME-circle arc fragments
      if (minHits > 0)
         praTask->SetMinNumHits(minHits);
      praTask->SetChargeFromCenter(true); // robust charge sign (angular sweep about the circle centre)
   }
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   // ─── UKF: pi+/pi- hypotheses (same PUMA calibration as run_digi_attpc.C) ───
   {
      const double e_C = 1.602176634e-19;
      const double u_MeV = 931.49410372;
      // hypothesis list: 2 (single species) or 4 (species="both": pi+/pi-/K+/K-)
      struct Hyp { std::string name; int sign; double mass; };
      std::vector<Hyp> hyps;
      if (isBoth) {
         hyps = {{"pi+", +1, 139.57039}, {"pi-", -1, 139.57039},
                 {"K+", +1, 493.677},    {"K-", -1, 493.677}};
      } else {
         hyps = {{nm + "+", +1, mMeV}, {nm + "-", -1, mMeV}};
      }

      auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
      for (const auto &hyp : hyps) {
         auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
         eloss->SetProjectile(1, 1, hyp.mass / u_MeV);
         std::vector<std::tuple<int, int, int>> mat;
         mat.push_back(std::make_tuple(18, 40, 9)); // Ar
         mat.push_back(std::make_tuple(6, 12, 1));  // C (from CH4)
         mat.push_back(std::make_tuple(1, 1, 4));   // H
         eloss->SetMaterial(mat);

         auto ukf = std::make_unique<EventFit::AtFitterUKF>(hyp.sign * e_C, hyp.mass, std::move(eloss));
         ukf->SetBField(ROOT::Math::XYZVector(0, 0, 4.0));
         ukf->SetUKFParameters(1e-3, 2.0, 0.0);
         ukf->SetMeasurementSigma(0.5);
         ukf->SetMomentumSigmaFrac(0.1);
         ukf->SetEnableEnergyStraggling(true);
         ukf->SetMinClusters(2);
         ukf->SetNIterations(3);
         ukf->SetZPadPlane(0.0);
         ukf->SetBackExtrapMaxPath(250.0);
         ukf->SetUseHelixBackExtrap(true);
         ukf->SetMinClusterSpacing(1.0);
         ukf->SetTargetClusters(clusterTarget);
         ukf->SetAdaptiveDistBounds(4.0, 14.0);
         ukf->SetUseClusterDirSeed(true);
         ukf->SetUseArcWalk(true);
         if (ringClustering) { // resistive ring centroiding + near-straight guard
            ukf->SetUseRingClustering(true, nPads);
            ukf->SetMaxSeedRadius(1500.0);
         }
         ukf->SetVertexZBias(useMFimpulse ? 0.0 : 8.6); // no shaping bias to correct in impulse mode
         // Deterministic vertex material correction (Cu trap + Al cryostat) — same knob
         // as genfit, so BOTH fitters report the at-vertex (not in-gas) momentum.
         if (matCorrRefDp > 0)
            ukf->SetVertexMaterialCorrection(375.0, matCorrRefDp);
         multi->AddHypothesis(std::move(ukf), hyp.name, hyp.sign);
      }
      auto *ukfTask = new AtFitterTask(std::move(multi));
      ukfTask->SetInputBranch("AtPatternEvent");
      ukfTask->SetOutputBranch("AtTrackingEventUKF");
      ukfTask->SetFitMetadataBranch("AtFitMetadataUKF");
      ukfTask->SetPersistence(kTRUE);
      fRun->AddTask(ukfTask);
   }

   // ─── genfit: pion, native Bethe-Bloch, per-track charge from PRA curvature ───
   {
      const double u_MeV = 931.49410372;
      // pdgPos (K+ / pi+); sign flipped per track from PRA charge.
      // matEffects ON: Bethe-Bloch forces effects on regardless of this flag, but
      // pass kFALSE so the intent (use material effects) is explicit.
      auto genfitter = std::make_unique<EventFit::AtGenfitter>(/*B*/ 4.0, pdgPos, mMeV / u_MeV, /*Z*/ 1,
                                                               /*eLossFile*/ "", /*noMatEffects*/ kFALSE, 2, 5);
      genfitter->SetEnergyLossBetheBloch(true);   // native Bethe-Bloch (pion-correct)
      genfitter->SetUseTrackChargeSign(true);     // pi+/pi- per PRA curvature
      genfitter->SetChargeSignFlip(genfitChargeFlip);
      genfitter->SetZPadPlane(0.0);               // match the PUMA UKF frame
      genfitter->SetZDriftSign(+1.0);             // PUMA sim hits are already lab-positive z
                                                  // (in the gas volume z=[0,300] mm); do NOT flip
      // Ring centroids are sub-pad (charge-weighted over ~several azimuthal pads),
      // so the measurement error is well below one pad; raw arc-walk clusters are
      // ~1 pad. Tell genfit the truth so it doesn't over-inflate the measurement cov.
      genfitter->SetMeasSigma(ringClustering ? 0.3 : 1.0);
      // PRA leaves ~2 clusters/track on PUMA's fine pad plane; re-cluster raw hits
      // into ~8 (same as the UKF) so genfit has enough points for a curvature fit.
      genfitter->SetReclusterArcWalk(true, clusterTarget, 3, 10);
      if (ringClustering) { // resistive ring centroiding + near-straight guard (no RK crash)
         genfitter->SetUseRingClustering(true, nPads);
         genfitter->SetMaxSeedRadius(1500.0); // skip near-straight tracks (p>1.8 GeV) before the fit
      }
      // Back-extrapolate the fitted state to the beam-axis POCA so the genfit vertex
      // is directly comparable to the UKF's (both report the at-vertex state).
      genfitter->SetBackExtrapPOCA(true);
      // Optionally keep material effects on in that extrap so the at-vertex momentum
      // recovers the ~20 MeV/c the pion lost in the Cu trap + Al cryostat before the gas.
      genfitter->SetBackExtrapMaterial(backExtrapMat);
      // Deterministic material correction (recover vertex p): calibrated for PUMA at
      // 375 MeV/c where the Cu/Al loss is ~18.2 MeV/c, scaled by Bethe-Bloch. Opt-in.
      if (matCorrRefDp > 0)
         genfitter->SetVertexMaterialCorrection(375.0, matCorrRefDp);
      // Branch-8 pions are in-plane (pz~0) so drift-z is degenerate for ordering;
      // pick the vertex end by beam-axis distance instead (fixes the vertex dr-tail).
      genfitter->SetVertexByXYRadius(true);
      genfitter->SetVertexZBias(useMFimpulse ? 0.0 : 8.6); // no shaping bias to correct in impulse mode
      genfitter->SetMinClusters(3);
      genfitter->SetVertexAxisMaxDist(150.0);
      auto *genfitTask = new AtFitterTask(std::move(genfitter));
      genfitTask->SetInputBranch("AtPatternEvent");
      genfitTask->SetOutputBranch("AtTrackingEventGenfit");
      genfitTask->SetFitMetadataBranch("AtFitMetadataGenfit");
      genfitTask->SetPersistence(kTRUE);
      fRun->AddTask(genfitTask);
   }

   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\n\033[1;32mrun_digi_ukf_genfit_test8 done.\033[0m\n"
             << "  out: " << outputFile << "\n"
             << "  branches: AtTrackingEventUKF, AtTrackingEventGenfit (+ AtPatternEvent, AtEvent)\n"
             << "  Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
