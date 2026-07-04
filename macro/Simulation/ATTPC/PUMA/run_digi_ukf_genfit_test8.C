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
                               TString psaType = "mfimpulse", TString species = "pi")
{
   const bool isK = (species == "K" || species == "kaon");
   const double mMeV = isK ? 493.677 : 139.57039; // kaon / pion mass
   const int pdgPos = isK ? 321 : 211;            // K+ / pi+
   const std::string nm = isK ? "K" : "pi";

   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi_both8.root";
   TString paramFile = "ATTPC.PUMA_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";
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
   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   mapping->GeneratePadPlane();

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   // PSA: "max" = AtPSAMax (peak; carries the +8.6 mm shaping-delay z bias, corrected
   // downstream via SetVertexZBias). "mfimpulse" = AtPSAMultiFit fitting the SAME GET
   // response the digi used, with z taken from the fitted IMPULSE time t0 -> removes the
   // shaping-delay z offset AT THE SOURCE (no SetVertexZBias needed).
   const bool useMFimpulse = (psaType == "mfimpulse");
   AtPSAtask *psaTask = nullptr;
   if (useMFimpulse) {
      auto mf = std::make_unique<AtPSAMultiFit>();
      mf->SetThreshold(0);
      mf->SetPeakingTime(0.5);       // 500 ns, matches PUMA AtPulse shaping (AtNominalResponse)
      mf->SetUseImpulseTime(true);   // z from fitted t0, not the response peak
      psaTask = new AtPSAtask(std::move(mf));
   } else {
      auto psa = std::make_unique<AtPSAMax>();
      psa->SetThreshold(0);
      psaTask = new AtPSAtask(std::move(psa));
   }
   psaTask->SetPersistence(kTRUE);

   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetTcluster(tCluster);
   praTask->SetTCUseSelectAndMerge(false);
   // Robust charge sign: angular sweep about the fitted circle centre (huge lever arm)
   // instead of the 3-point cross product, which is noise-dominated on PUMA's shallow
   // 312 mm-radius arcs. Lifts charge accuracy 83% -> ~99.6% for BOTH fitters.
   praTask->SetChargeFromCenter(true);
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   // ─── UKF: pi+/pi- hypotheses (same PUMA calibration as run_digi_attpc.C) ───
   {
      const double e_C = 1.602176634e-19;
      const double u_MeV = 931.49410372;
      const int signs[2] = {+1, -1};
      const std::string names[2] = {nm + "+", nm + "-"};

      auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
      for (int i = 0; i < 2; ++i) {
         auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
         eloss->SetProjectile(1, 1, mMeV / u_MeV);
         std::vector<std::tuple<int, int, int>> mat;
         mat.push_back(std::make_tuple(18, 40, 9)); // Ar
         mat.push_back(std::make_tuple(6, 12, 1));  // C (from CH4)
         mat.push_back(std::make_tuple(1, 1, 4));   // H
         eloss->SetMaterial(mat);

         auto ukf = std::make_unique<EventFit::AtFitterUKF>(signs[i] * e_C, mMeV, std::move(eloss));
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
         ukf->SetTargetClusters(8);
         ukf->SetAdaptiveDistBounds(4.0, 14.0);
         ukf->SetUseClusterDirSeed(true);
         ukf->SetUseArcWalk(true);
         ukf->SetVertexZBias(useMFimpulse ? 0.0 : 8.6); // no shaping bias to correct in impulse mode
         multi->AddHypothesis(std::move(ukf), names[i], signs[i]);
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
      genfitter->SetMeasSigma(1.0);               // near the digi resolution floor
      // PRA leaves ~2 clusters/track on PUMA's fine pad plane; re-cluster raw hits
      // into ~8 (same as the UKF) so genfit has enough points for a curvature fit.
      genfitter->SetReclusterArcWalk(true, 8, 3, 10);
      // Back-extrapolate the fitted state to the beam-axis POCA so the genfit vertex
      // is directly comparable to the UKF's (both report the at-vertex state).
      genfitter->SetBackExtrapPOCA(true);
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
