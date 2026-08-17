/// @file fitGenfitter_a1975.C
/// @brief Fit a pre-reconstructed a1975 file with the CLEAN EventFit::AtGenfitter
/// (modern AtFitter base + AtFitterTask, geometric closest-to-axis seed, no merging).
/// Reads <run>_reco.root (AtPatternEvent), writes <run>_genfitter<suffix>.root with
/// AtTrackingEvent + AtPIDEvent. Sibling of fitGenfit_a1975.C (the legacy fallback).
///
///   root -b -q 'pipeline/fitGenfitter_a1975.C("run_0106", 500, "/mnt/f/a1975/reco/", "_test", "/tmp/")'

void fitGenfitter_a1975(TString fileName = "run_0106", Long64_t nEvents = -1, TString ioDir = "/mnt/f/a1975/reco/",
                        TString outSuffix = "", TString outDir = "", Double_t bField = -2.85, Int_t minIter = 2,
                        Int_t maxIter = 5, TString elossName = "deuteron_H2_catima.txt",
                        TString pidGate = "pid/deuteron_band.json", Double_t measSigma = 4.0,
                        Bool_t matEffects = kFALSE, Int_t pdg = 1000010020, Double_t massAmu = 2.0135532,
                        Int_t Z = 1, Bool_t mergeContinuity = kFALSE,
                        // Geometry is a PARAMETER, not a constant: the a1975 gas is H2 at 300
                        // torr (3.308e-5 g/cm3) and this defaulted to ATTPC_H1bar (8.27e-5),
                        // i.e. 2.5x too much material. Whether that matters depends on
                        // matEffects -- with it off genfit applies no material effects and the
                        // density cannot reach the fit. Left at the historical default so an
                        // existing call reproduces the existing production; pass the correct one
                        // explicitly to test or to redo.
                        TString geoName = "ATTPC_H1bar", Bool_t catimaMSC = kFALSE,
                        Bool_t catimaStraggling = kFALSE, Bool_t matFallback = kTRUE,
                        Double_t elossDensityGCm3 = 0,
                        // getFittedState() is the FIRST MEASUREMENT POINT, not the vertex: the
                        // ejectile already crossed unmeasured gas before its first cluster and
                        // that energy is missing from the fit. backExtrap transports the state
                        // back to the beam axis, which recovers it -- but ONLY while matEffects
                        // is on, since without material the transport is geometric and |p| comes
                        // back unchanged. It fills GetKinematicsXtr(); GetKinematics() keeps the
                        // raw fit either way.
                        //
                        // There is deliberately NO manual dE/dx knob here. AtGenfitter's
                        // SetManualELoss applies its own loss over the gap unconditionally, so
                        // combining it with matEffects counts the gap TWICE (measured on (d,t):
                        // +9.5% -> +4.65% on the low branch once the manual term was dropped).
                        // Last param, so every existing positional caller is unaffected.
                        Bool_t backExtrap = kFALSE,
                        // Take dE/dx from CATIMA instead of the ASCII curve loaded above. The
                        // table is consulted ONLY below beta*gamma = 0.05 (KE 2.3 MeV for a
                        // deuteron, 3.5 for a triton) and is multiplied by one GLOBAL density,
                        // while Bethe-Bloch above that uses the per-step material -- so the two
                        // branches disagree about what they are traversing. CATIMA uses the step's
                        // own material and needs no table at all.
                        // catimaELossFull additionally replaces Bethe-Bloch; validate separately.
                        // Expect little effect on (p,d): these deuterons are ~20 MeV and only
                        // touch the table region near the endpoint. It is (d,t), whose low branch
                        // is 0.8-6 MeV against a 3.5 MeV threshold, that sits inside it.
                        Bool_t catimaELoss = kFALSE, Bool_t catimaELossFull = kFALSE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/" + geoName + "_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter" + outSuffix + ".root";
   // empty elossName -> empty eloss file: rely on genfit's internal Bethe-Bloch (material effects)
   std::string elossFile = elossName.Length() ? ((std::string)dir.Data() + "/resources/energy_loss/" + elossName.Data()) : std::string("");

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfitter_a1975 (clean) ===\033[0m\n  in  : " << inputFile << "\n  out : "
             << outputFile << "\n  B=" << bField << " iter " << minIter << "-" << maxIter << "\n";

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

   // particle hypothesis from args (6th ctor arg is noMatEffects -> invert matEffects)
   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << "\n";
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, elossFile, !matEffects, minIter,
                                                         maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);

   // ---- CATIMA material model ---------------------------------------------------------------
   // Only reachable with matEffects = kTRUE: noiseCoulomb / noiseBetheBloch are the sole places
   // the backend is consulted and neither is called with material effects off.
   //
   // TURNING matEffects ON MAKES THE GEOMETRY LIVE. With it off the gas density cannot reach the
   // fit, which is why this macro's historical geoName default (ATTPC_H1bar, 8.27e-5 g/cm3) was
   // harmless. The a1975 gas is H2 at 300 torr: ATTPC_H300torr_RT is 3.308e-5, i.e. the default
   // is 2.5x too much material. Pass the right geometry whenever matEffects is on.
   if (catimaMSC || catimaStraggling) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catima flags set but matEffects is OFF -- they are inert.\033[0m\n";
      if (geoName.Contains("H1bar"))
         std::cout << "\033[1;31mWARNING: matEffects with geoName=" << geoName
                   << " -- that gas is 8.27e-5 g/cm3, 2.5x the a1975 H2 at 300 torr. "
                      "Use ATTPC_H300torr_RT.\033[0m\n";
      fitter->SetCatimaMaterial(catimaMSC, catimaStraggling);
      std::cout << "  \033[1;35mCATIMA material model: MSC " << (catimaMSC ? "ON" : "off")
                << ", straggling " << (catimaStraggling ? "ON" : "off") << "\033[0m\n";
   }
   // Hybrid dE/dx: the table serves beta*gamma < 0.05 where genfit would otherwise apply ZERO
   // stopping power. For a deuteron that threshold is KE = 2.3 MeV.
   if (elossDensityGCm3 > 0 && !elossFile.empty()) {
      fitter->SetELossHybrid(kTRUE, elossDensityGCm3);
      std::cout << "  \033[1;32mdE/dx TABLE in hybrid mode, density " << elossDensityGCm3
                << " g/cm3\033[0m\n";
   }
   // A throwing track must NOT be silently refitted without material and kept: that blends two
   // physics models into one spectrum. Off for any material-effects production.
   fitter->SetMatEffectsFallback(matFallback);
   // Back-extrapolation to the beam axis. Announced rather than silent, because whether it did
   // anything depends on matEffects: with material off it is a no-op on |p| and Xtr comes back
   // equal to the raw fit, which looks identical to having forgotten the flag.
   // CATIMA dE/dx. Announced only when ON, and it warns rather than silently doing nothing if
   // material effects are off -- with matEffects off no dE/dx model is consulted at all, so the
   // flag would be inert and the run would look like a valid arm of an A/B when it is not.
   fitter->SetCatimaELoss(catimaELoss, catimaELossFull);
   if (catimaELoss) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catimaELoss set but matEffects is OFF -- it is inert.\033[0m\n";
      else
         std::cout << "  \033[1;32mdE/dx from CATIMA"
                   << (catimaELossFull ? " over the FULL range (Bethe-Bloch replaced too)"
                                       : " below beta*gamma=0.05 (Bethe-Bloch kept above)")
                   << ", per-step material\033[0m\n";
   }

   fitter->SetBackExtrapToAxis(backExtrap);
   if (backExtrap) {
      std::cout << "  \033[1;32mback-extrapolation to the beam axis: ON"
                << (matEffects ? " (recovers the vertex-gap energy)"
                               : " -- but matEffects is OFF, so |p| will come back UNCHANGED")
                << "\033[0m\n";
   }
   if (mergeContinuity) {
      fitter->SetMergeContinuity(kTRUE); // merge PRA-split fragments of one track before fitting
      std::cout << "  continuity merging: ON\n";
   }
   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating will be disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "\n";
      }
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
}
