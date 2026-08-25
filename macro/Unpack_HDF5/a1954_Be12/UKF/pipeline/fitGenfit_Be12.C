/// @file fitGenfit_Be12.C
/// @brief Fit a pre-reconstructed a1954 12Be(p,p') file with EventFit::AtGenfitter
///        (updated GenFit fork: Yassid/GenFit, SRIM energy-loss tables).
///
/// Ported from a1975 D2_UKF/fitGenfitter_a1975_deuterium.C. Reads <run>_reco.root
/// (AtPatternEvent) and writes <run>_genfit<suffix>.root (AtTrackingEvent + AtPIDEvent).
/// Same input as fitUKF_Be12.C, so the two fitters can be compared track-by-track.
///
/// >>> REQUIRES the GENFIT-enabled build:  source build_genfit/config.sh  <<<
///     (the default build/ has GENFIT2 NOT FOUND, so AtGenfitter is a stub there.)
///
/// HANDEDNESS: experimental data -> Bz negative (as for the UKF). Default bField=-2.85;
///     verify on fitted vertex/KE and flip if fits diverge (bField is an arg).
///
///   root -b -q 'fitGenfit_Be12.C("run_0142", 500, "/home/yassid/a1954_Be12_reco/", "_test", "/tmp/")'
///   root -b -q 'fitGenfit_Be12.C("run_0142")'   // whole run, in-place output

void fitGenfit_Be12(TString fileName = "run_0142", Long64_t nEvents = -1,
                    TString ioDir = "/home/yassid/a1954_Be12_reco/", TString outSuffix = "", TString outDir = "",
                    Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5, TString pidGate = "",
                    Double_t measSigma = 4.0, Double_t thetaMinDeg = 10.0, Double_t thetaMaxDeg = 170.0,
                    // matEffects DEFAULT FLIPPED kFALSE -> kTRUE 2026-08-25: the a1975 working
                    // point adopts material effects everywhere, and they are only safe now that
                    // the geometry below is the real 300-torr gas. genfit_hdb_batch.sh calls this
                    // with defaults only, so the default is what that production gets.
                    Bool_t matEffects = kTRUE, Bool_t backwardSeedFix = kFALSE, TString particle = "proton",
                    TString geoName = "ATTPC_H300torr_RT",
                    // ---- material model (ported from fitGenfitter_a1975.C, 2026-08-25) -------
                    // matEffects ON is only safe with ALL FOUR co-requisites, each of which fails
                    // SILENTLY on its own (see a1975 working_point.sh):
                    //   correct geometry .... ATTPC_H300torr_RT above, now the default
                    //   matFallback FALSE ... a throwing track must NOT be silently refitted
                    //                         without material and kept; that blends two physics
                    //                         models into one spectrum
                    //   no manual dE/dx ..... AtGenfitter::SetManualELoss runs unconditionally, so
                    //                         pairing it with matEffects double-counts the vertex
                    //                         gap. There is deliberately no knob for it here.
                    //   a dE/dx source ...... genfit applies ZERO stopping power below
                    //                         beta*gamma = 0.05. For a PROTON that is KE 1.17 MeV
                    //                         (3.4 % of the a1954 (p,p') tracks); for a deuteron
                    //                         2.3 and a triton 3.5. Supply it with catimaELoss
                    //                         (per-step, per-material, no table) or with
                    //                         elossName + elossDensity (frozen ASCII table).
                    // ADOPTED on a1975 after the (d,t) A/B: catimaELoss = kTRUE, full range OFF --
                    // collapsed fits 2.51 % -> 1.08 % with energy scale and resolution unchanged.
                    Bool_t catimaELoss = kTRUE, Bool_t catimaELossFull = kFALSE,
                    Bool_t catimaMSC = kFALSE, Bool_t catimaStraggling = kFALSE,
                    Bool_t matFallback = kFALSE, TString elossName = "",
                    Double_t elossDensityGCm3 = 0,
                    // getFittedState() is the FIRST MEASUREMENT POINT, not the vertex: the ejectile
                    // crossed unmeasured gas before its first cluster. backExtrap transports the
                    // state back to the beam axis and recovers it -- but ONLY while matEffects is
                    // on; with material off the transport is geometric and |p| returns unchanged.
                    // Fills GetKinematicsXtr(); GetKinematics() keeps the raw fit either way.
                    Bool_t backExtrap = kFALSE)
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
   // a1954 ran H2 at 300 torr (CONFIRMED 2026-08-25) -> ATTPC_H300torr_RT, rho = 3.308e-5
   // g/cm3 (media.geo H_300torr_RT). The old default ATTPC_H600torr is 6.616e-5, i.e. TWICE
   // the material; it is inert while matEffects = kFALSE but goes live the moment it is on.
   TString geoManFile = dir + "/geometry/" + geoName + "_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954_Be12.par";
   // empty elossName -> empty eloss file: dE/dx comes from CATIMA (catimaELoss) or from
   // genfit's internal Bethe-Bloch. Resolved from resources/energy_loss/ as in a1975.
   std::string elossFile =
      elossName.Length() ? ((std::string)dir.Data() + "/resources/energy_loss/" + elossName.Data()) : std::string("");
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfit" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found. Run unpackReco_Be12.C first.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfit_Be12 (12Be(p,p'), PROTON hyp) ===\033[0m\n"
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
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, elossFile, !matEffects, minIter,
                                                         maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);

   // ---- material model ------------------------------------------------------------------------
   // Everything below is reachable ONLY with matEffects = kTRUE: noiseCoulomb / noiseBetheBloch are
   // the sole places the backend is consulted and neither is called with material effects off. Each
   // flag therefore WARNS instead of silently doing nothing, so a run cannot look like a valid arm
   // of an A/B when it is inert.
   if (matEffects && geoName.Contains("600torr"))
      std::cout << "\033[1;31mWARNING: matEffects with geoName=" << geoName
                << " -- that gas is 6.616e-5 g/cm3, TWICE the a1954 H2 at 300 torr. "
                   "Use ATTPC_H300torr_RT.\033[0m\n";
   if (catimaMSC || catimaStraggling) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catima material flags set but matEffects is OFF -- inert.\033[0m\n";
      fitter->SetCatimaMaterial(catimaMSC, catimaStraggling);
      std::cout << "  \033[1;35mCATIMA material model: MSC " << (catimaMSC ? "ON" : "off") << ", straggling "
                << (catimaStraggling ? "ON" : "off") << "\033[0m\n";
   }
   // Hybrid ASCII dE/dx: the table serves beta*gamma < 0.05, where genfit applies ZERO stopping.
   if (elossDensityGCm3 > 0 && !elossFile.empty()) {
      fitter->SetELossHybrid(kTRUE, elossDensityGCm3);
      std::cout << "  \033[1;32mdE/dx TABLE in hybrid mode, density " << elossDensityGCm3 << " g/cm3\033[0m\n";
   }
   fitter->SetMatEffectsFallback(matFallback);
   if (matEffects && matFallback)
      std::cout << "\033[1;31mWARNING: matFallback is ON -- throwing tracks get refitted WITHOUT "
                   "material and kept, mixing two physics models in one spectrum.\033[0m\n";
   fitter->SetCatimaELoss(catimaELoss, catimaELossFull);
   if (catimaELoss) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catimaELoss set but matEffects is OFF -- inert.\033[0m\n";
      else
         std::cout << "  \033[1;32mdE/dx from CATIMA"
                   << (catimaELossFull ? " over the FULL range (Bethe-Bloch replaced too)"
                                       : " below beta*gamma=0.05 (Bethe-Bloch kept above)")
                   << ", per-step material\033[0m\n";
   } else if (matEffects && elossFile.empty()) {
      std::cout << "\033[1;31mWARNING: matEffects ON with NO dE/dx source -- genfit applies zero "
                   "stopping power below beta*gamma=0.05 (KE 1.17 MeV for a proton).\033[0m\n";
   }
   fitter->SetBackExtrapToAxis(backExtrap);
   if (backExtrap)
      std::cout << "  \033[1;32mback-extrapolation to the beam axis: ON"
                << (matEffects ? " (recovers the vertex-gap energy)"
                               : " -- but matEffects is OFF, so |p| comes back UNCHANGED")
                << "\033[0m\n";
   // --------------------------------------------------------------------------------------------

   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);
   fitter->SetBackwardSeedFix(backwardSeedFix);
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
}
