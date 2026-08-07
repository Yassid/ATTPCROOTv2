/// @file fitGenfit_C14.C
/// @brief Fit a pre-reconstructed a1954 14C(p,p') file with EventFit::AtGenfitter
///        (updated GenFit fork: Yassid/GenFit, SRIM energy-loss tables).
///
/// Ported from a1975 D2_UKF/fitGenfitter_a1975_deuterium.C. Reads <run>_reco.root
/// (AtPatternEvent) and writes <run>_genfit<suffix>.root (AtTrackingEvent + AtPIDEvent).
/// Same input as fitUKF_C14.C, so the two fitters can be compared track-by-track.
///
/// >>> REQUIRES the GENFIT-enabled build:  source build_genfit/config.sh  <<<
///     (the default build/ has GENFIT2 NOT FOUND, so AtGenfitter is a stub there.)
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
                    TString geoName = "ATTPC_H300torr", Bool_t matFallback = kTRUE,
                    Bool_t backExtrap = kFALSE, Double_t manualElossDensity = 0.0, Int_t matA = 1)
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
   if (manualElossDensity > 0) {
      fitter->SetManualELoss(manualElossDensity, matA);
      std::cout << "  Manual CATIMA eloss over the vertex gap: rho=" << manualElossDensity
                << " g/cm3, target A=" << matA << std::endl;
   }
   fitter->SetBackwardSeedFix(backwardSeedFix);
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
}
