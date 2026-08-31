/// @file fitGenfit_C14.C
/// @brief Fit a pre-reconstructed a1954 14C(p,p') file with EventFit::AtGenfitter
///        (updated GenFit fork: Yassid/GenFit, SRIM energy-loss tables).
///
/// Ported from a1975 D2_UKF/fitGenfitter_a1975_deuterium.C. Reads <run>_reco.root
/// (AtPatternEvent) and writes <run>_genfit<suffix>.root (AtTrackingEvent + AtPIDEvent).
/// Same input as fitUKF_C14.C, so the two fitters can be compared track-by-track.
///
/// Build: `source build/config.sh`. The two-build split this file used to require
/// (`build_genfit/` for GENFIT, `build/` for everything else) is gone -- there is now a single
/// build and it is GENFIT-enabled: build/CMakeCache.txt resolves GENFIT2_LIBRARY to
/// ~/fair_install/GenFit/lib/libgenfit2.so and libAtReconstruction.so exports a real
/// EventFit::AtGenfitter, not a stub. pp/refit_genfit_xtr.sh already sources build/config.sh.
///
/// HANDEDNESS: experimental data -> Bz negative (as for the UKF). Default bField=-2.85;
///     verify on fitted vertex/KE and flip if fits diverge (bField is an arg).
///
///   root -b -q 'fitGenfit_C14.C("run_0055", 500, "/home/yassid/a1954_C14_reco/", "_test", "/tmp/")'
///   root -b -q 'fitGenfit_C14.C("run_0055")'   // whole run, in-place output

void fitGenfit_C14(TString fileName = "run_0055", Long64_t nEvents = -1,
                    TString ioDir = "/home/yassid/a1954_C14_reco/", TString outSuffix = "", TString outDir = "",
                    Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5, TString pidGate = "",
                    Double_t measSigma = 4.0, Double_t thetaMinDeg = 10.0, Double_t thetaMaxDeg = 170.0,
                    Bool_t matEffects = kFALSE, Bool_t backwardSeedFix = kFALSE, TString particle = "proton",
                    // ATTPC_H300torr is H2 at 300 torr and 0 C (3.553e-5). The a1954 gas is at ROOM
                    // TEMPERATURE: ATTPC_H300torr_RT, 3.308e-5 -- 7.4 % less material, and the value
                    // the adopted gf_xtr production already used for its manual eloss. Inert while
                    // matEffects is off; LIVE the moment it is on, which is why the default moves.
                    TString geoName = "ATTPC_H300torr_RT", Bool_t matFallback = kTRUE,
                    Bool_t backExtrap = kFALSE, Double_t manualElossDensity = 0.0, Int_t matA = 1,
                    // ---- CATIMA (added 2026-08-25) ------------------------------------------
                    // Adopted setting on a1975 after the (d,t) A/B: catimaELoss ON, full range OFF
                    // -- collapsed fits 2.51 % -> 1.08 % with energy scale and resolution unchanged.
                    // Only reachable with matEffects = kTRUE; each flag warns rather than silently
                    // doing nothing, so an inert run cannot pass for a valid arm of a comparison.
                    Bool_t catimaELoss = kFALSE, Bool_t catimaELossFull = kFALSE,
                    Bool_t catimaMSC = kFALSE, Bool_t catimaStraggling = kFALSE,
                    // ---- RANGE CONSTRAINT (exposed 2026-08-29) ------------------------------
                    // For a track that STOPS in the gas, the path length measures the energy far
                    // better than the curvature does. That matters for backward (d,p) protons:
                    // they spiral several turns, shed a large part of their energy on the way, and
                    // the pattern circle then returns the AVERAGE radius rather than the one at the
                    // vertex -- 0.895 of truth, which is 0.80 in energy and exactly the -20 % bias
                    // those tracks show. AtGenfitter has had the machinery (containment + a Bragg
                    // dQ/dx test + straggling errors) all along; it was simply not reachable from
                    // here. rangeDensity is g/cm3 and rangeMatA is 1 for H2, 2 for D2.
                    Bool_t rangeConstraint = kFALSE, Double_t rangeDensity = 0.0, Int_t rangeMatA = 2,
                    // ---- FIRST-ARC SEED (exposed 2026-08-29) --------------------------------
                    // The default seed takes the momentum from GetGeoRadius(), a single circle
                    // fitted to the WHOLE track. For a backward proton that spirals to a stop that
                    // circle is the average of a tightening spiral -- measured at 0.895 of the
                    // radius at the vertex, i.e. 0.80 in energy, which is the -20 % bias those
                    // tracks carry. AtSpyralPID instead fits the circle to the FIRST ARC and
                    // regresses rho against z, which is the radius the analysis actually wants.
                    // Building the estimator needs no PID gate despite what the header says.
                    Bool_t seedFromSpyral = kFALSE, Bool_t useClusterOrder = kFALSE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   // ejectile hypothesis: proton for (p,p'), deuteron for (p,d)
   Int_t pdg = 2212;
   Double_t massAmu = 1.00782503207;
   Int_t Z = 1;
   if (particle == "deuteron") {
      pdg = 1000010020;
      massAmu = 2.01355321275;
      Z = 1;
   } else if (particle == "triton") {
      pdg = 1000010030;
      massAmu = 3.01550071632;
      Z = 1;
   }
   std::cout << "particle hypothesis: " << particle << "\n";

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   // a1954 ran H2 at 300 torr -> ATTPC_H300torr (rho = 3.553e-5 g/cm3). The first production
   // used ATTPC_H600torr (6.616e-5), ~1.9x too dense. With matEffects = kFALSE (the default)
   // the geometry only drives navigation, so that mistake was a no-op for genfit -- verified:
   // the 600- and 300-torr genfit caches are identical to the last digit. It matters only if
   // material effects are switched on.
   //
   // matEffects = kFALSE is the DEFAULT ON PURPOSE. Turning it on does fix the absolute scale
   // (g.s. +0.141 -> +0.019, agreeing with UKF to 0.071 MeV) but on a1954 14C it also costs
   // 32% of the tracks, degrades FWHM 0.364 -> 0.487 and triples the angle drift, and it is
   // ~20-45x slower because unstable tracks (ill-conditioned covariance / Cholesky failures)
   // are silently REFIT WITHOUT material effects by AtGenfitter, mixing two populations into
   // one output with no flag. Fix that instability before making it the default.
   TString geoManFile = dir + "/geometry/" + geoName + "_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954_C14.par";
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfit" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found. Run unpackReco_C14.C first.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfit_C14 (14C(p,p'), PROTON hyp) ===\033[0m\n"
             << "  in   : " << inputFile << "\n  out  : " << outputFile << "\n  B=" << bField << "  iter " << minIter
             << "-" << maxIter << "  theta[" << thetaMinDeg << "," << thetaMaxDeg << "]\n";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   // genfit material/propagation needs gGeoManager loaded before the fit
   if (gROOT->FindObject("FAIRGeom") == nullptr) {
      TFile *gf = TFile::Open(geoManFile);
      gf->Get("FAIRGeom");
   }

   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << " (proton)\n";
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, /*elossFile*/ "", !matEffects,
                                                         minIter, maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);
   // Vertex-gap energy recovery, the a1975 D2 pattern (fitGenfitter_a1975_deuterium.C:93-102).
   // genfit's getFittedState() is the FIRST MEASUREMENT POINT: the gas between the reaction
   // vertex and the first cluster is unmeasured, so its energy loss is missing from |p|. The
   // extrapolation supplies the gap length geometrically and CATIMA supplies the dE/dx over it,
   // which works with matEffects OFF -- no need to pay the ~60 % good-fit rate that turning
   // genfit material effects on costs. Result lands in GetKinematicsXtr(), NOT GetKinematics().
   // matA = 1 for the H2 target here (2 would be deuterium).
   fitter->SetBackExtrapToAxis(backExtrap);
   if (backExtrap)
      std::cout << "  Back-extrapolation to beam axis: ON" << std::endl;
   // DOUBLE-COUNT GUARD. SetManualELoss applies its own loss over the vertex gap
   // UNCONDITIONALLY. With matEffects ON, extrapolateToLine already integrates that gap, so the
   // two together count it TWICE -- measured on a1975 (d,t) as +9.5 % falling to +4.65 % once the
   // manual term was dropped. The manual term is the matFX-OFF substitute for material effects,
   // not a companion to them.
   if (manualElossDensity > 0 && matEffects) {
      std::cout << "\033[1;31mREFUSING manual eloss: matEffects is ON, so extrapolateToLine already "
                   "integrates the vertex gap and this would count it twice. Ignored.\033[0m\n";
   } else if (manualElossDensity > 0) {
      fitter->SetManualELoss(manualElossDensity, matA);
      std::cout << "  Manual CATIMA eloss over the vertex gap: rho=" << manualElossDensity
                << " g/cm3, target A=" << matA << std::endl;
   }
   // ---- CATIMA material model ---------------------------------------------------------------
   if (matEffects && geoName.Contains("H300torr") && !geoName.Contains("_RT"))
      std::cout << "\033[1;31mWARNING: matEffects with " << geoName << " = 3.553e-5 g/cm3 (0 C). "
                   "The a1954 gas is 3.308e-5 at room temperature: use ATTPC_H300torr_RT.\033[0m\n";
   if (catimaMSC || catimaStraggling) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catima material flags set but matEffects is OFF -- inert.\033[0m\n";
      fitter->SetCatimaMaterial(catimaMSC, catimaStraggling);
      std::cout << "  \033[1;35mCATIMA material model: MSC " << (catimaMSC ? "ON" : "off")
                << ", straggling " << (catimaStraggling ? "ON" : "off") << "\033[0m\n";
   }
   fitter->SetCatimaELoss(catimaELoss, catimaELossFull);
   if (catimaELoss) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catimaELoss set but matEffects is OFF -- inert.\033[0m\n";
      else
         std::cout << "  \033[1;32mdE/dx from CATIMA"
                   << (catimaELossFull ? " over the FULL range (Bethe-Bloch replaced too)"
                                       : " below beta*gamma=0.05 (Bethe-Bloch kept above)")
                   << ", per-step material\033[0m\n";
   } else if (matEffects) {
      std::cout << "\033[1;31mWARNING: matEffects ON with NO dE/dx source -- genfit applies zero "
                   "stopping power below beta*gamma=0.05 (KE 1.17 MeV for a proton).\033[0m\n";
   }
   // A non-owning observer: `fitter` is moved into the task below, so anything read afterwards has
   // to go through this or it dereferences a moved-from unique_ptr (it crashed the first time).
   auto *fitterObs = fitter.get();
   if (seedFromSpyral) {
      fitter->SetSeedFromSpyral(kTRUE);
      std::cout << "  SEED FROM SPYRAL: circle on the first arc, polar from a rho-vs-z regression\n";
   }
   fitter->SetBackwardSeedFix(backwardSeedFix);
   // keep the PRA cluster order (pair with arc-walk) instead of re-sorting by drift z
   fitter->SetUseClusterOrder(useClusterOrder);
   if (rangeConstraint) {
      if (!(rangeDensity > 0))
         std::cout << "\033[1;31mWARNING: rangeConstraint on with density " << rangeDensity
                   << " -- the range<->energy table needs a real gas density, constraint will be inert.\033[0m\n";
      fitter->SetRangeConstraint(kTRUE, rangeDensity, rangeMatA);
      std::cout << "  RANGE CONSTRAINT on: rho = " << rangeDensity << " g/cm3, matA = " << rangeMatA << "\n";
   }
   // When a material-effects fit throws, the fitter retries it with material effects OFF.
   // Those tracks come from a different model and are marked in AtFitTrackMetadata
   // (GetMatEffects()==false, GetMatEffectsFallback()==true) so they can be filtered out.
   // Pass matFallback=kFALSE to drop them instead, giving a clean matFX-only sample.
   fitter->SetMatEffectsFallback(matFallback);
   if (matEffects)
      std::cout << "  matFX fallback: " << (matFallback ? "on (retries are flagged)" : "OFF (failures dropped)")
                << "\n";
   if (pidGate.Length() && !gSystem->AccessPathName(pidGate.Data())) {
      fitter->SetPIDGate(pidGate.Data());
      std::cout << "  PID gate: " << pidGate << "\n";
   } else {
      std::cout << "  PID gate: NONE (fit all tracks)\n";
   }

   AtFitterTask *fitterTask = new AtFitterTask(std::move(fitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   run->AddTask(pidTask);
   run->AddTask(fitterTask);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   t.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << t.RealTime() << " s)\n";
   // Say what the range constraint actually did. Without this the constraint can be enabled, run,
   // and change nothing, with no way to tell whether it was rejected by containment, by the Bragg
   // test, or never reached at all -- which is exactly what happened the first time it was tried.
   if (rangeConstraint && fitterObs)
      std::cout << "  range constraint: contained " << fitterObs->GetNRangeContained() << ", applied "
                << fitterObs->GetNRangeApplied() << ", failed the Bragg test " << fitterObs->GetNRangeBraggFail()
                << "\n";
}
