/// @file fitGenfitter_a1975_deuterium.C
/// @brief Fit a pre-reconstructed a1975 D2-target file (16C(d,p)17C) with the clean
/// EventFit::AtGenfitter, using the PROTON hypothesis (the detected ejectile in
/// (d,p) inverse kinematics). Extended to BACKWARD tracks: the proton goes largely
/// backward in the lab, so the theta acceptance window is 10-170 deg AND the backward
/// seed-reversal fix is ON by default (backwardSeedFix=kTRUE -> AtGenfitter::
/// SetBackwardSeedFix). Without it ~52% of backward (theta_geo>=90) protons get seeded
/// toward +z_lab and genfit reflects them into the forward hemisphere with wrong
/// theta/KE; the fix seeds them from the far (highest-z_lab) end. Backward yield
/// recovers ~3x (theta_lab>90 candidate fraction 5.6% -> 17.5% on run_0016). The
/// AtGenfitter framework default stays OFF, so forward (p,d)/(p,p) pipelines are unchanged.
///
/// Reads reco_d2/<run>_reco.root (AtPatternEvent), writes <run>_genfitter_p<suffix>.root
/// with AtTrackingEvent + AtPIDEvent. Sibling of UKF/pipeline/fitGenfitter_a1975.C
/// (which fits the DEUTERON hypothesis for the proton-target (p,d) channel).
///
///   root -b -q 'fitGenfitter_a1975_deuterium.C("run_0016", 500, "/mnt/f/a1975/reco_d2/", "_test", "/tmp/")'
///
/// B-sign note: a1975 is experimental data (z_lab = ZPadPlane - z_digi handedness).
/// The proton-target genfit (p,d) reference used B=+2.85; the deuterium unpackNFit
/// macro used -2.85. The correct sign for the D2 backward protons must be verified
/// on the fitted vertex/KE -- bField is an arg so both can be tried quickly.

/// Species is parametrised so the same macro serves every D2 ejectile; the defaults are the
/// PROTON hypothesis of 16C(d,p)17C, so existing calls are unchanged.
///   (d,p) 17C : defaults                          -> <run>_genfitter_p<suffix>.root
///   (d,t) 15C : pdg 1000010030, 3.01550072, Z 1, speciesTag "t"
/// The eloss file is empty either way (genfit's internal model), so no per-species table is needed.
void fitGenfitter_a1975_deuterium(TString fileName = "run_0016", Long64_t nEvents = -1,
                                  TString ioDir = "/mnt/f/a1975/reco_d2/", TString outSuffix = "",
                                  TString outDir = "", Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5,
                                  TString pidGate = "", Double_t measSigma = 4.0, Double_t thetaMinDeg = 10.0,
                                  Double_t thetaMaxDeg = 170.0, Bool_t matEffects = kFALSE,
                                  Bool_t backwardSeedFix = kTRUE, Int_t pdg = 2212,
                                  Double_t massAmu = 1.00782503207, Int_t Z = 1,
                                  TString speciesTag = "p", TString recoSuffix = "_reco",
                                  TString geoName = "ATTPC_H1bar_geomanager.root",
                                  Bool_t backExtrap = kFALSE, Double_t manualElossDensity = 0, Int_t matA = 2)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   // The geometry supplies genfit's MATERIAL, so it has to be the real target gas whenever
   // matEffects is on. The default here is the historical H1bar, whose H_1bar medium is
   // ordinary hydrogen and has nothing to do with a D2 target; pass
   // ATTPC_D300torr_v2_geomanager.root (medium TargetD2_300: A 2.0141, Z 1, 6.564e-5 g/cm3)
   // for these runs. With matEffects off the geometry only drives navigation and the choice
   // is harmless, which is why this went unnoticed.
   TString geoManFile = dir + "/geometry/" + geoName;
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par"; // D2 gas
   TString inputFile = ioDir + fileName + recoSuffix + ".root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter_" + speciesTag + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfitter_a1975_deuterium (D2 target, ejectile '" << speciesTag
             << "' pdg=" << pdg << " m=" << massAmu << " amu Z=" << Z << ") ===\033[0m\n"
             << "  in   : " << inputFile << "\n  out  : " << outputFile << "\n  B=" << bField << "  iter "
             << minIter << "-" << maxIter << "  theta[" << thetaMinDeg << "," << thetaMaxDeg << "]\n";

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

   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << " ('" << speciesTag << "')\n";
   // empty eloss file -> genfit internal; noMatEffects = !matEffects
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, /*elossFile*/ "", !matEffects,
                                                         minIter, maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg); // keep BACKWARD protons (default clips at 170)
   fitter->SetBackwardSeedFix(backwardSeedFix);      // prototype: fix backward-track seed reversal
   fitter->SetBackExtrapToAxis(backExtrap);
   if (backExtrap)
      std::cout << "  Back-extrapolation to beam axis: ON\n";
   // Hand-applied dE/dx over the vertex gap, for the matEffects-OFF configuration: genfit stays
   // fast and keeps its 85% good-fit rate, CATIMA supplies the energy the fit never saw.
   if (manualElossDensity > 0) {
      fitter->SetManualELoss(manualElossDensity, matA);
      std::cout << "  Manual CATIMA eloss over the vertex gap: rho=" << manualElossDensity
                << " g/cm3, target A=" << matA << "\n";
   }
   if (backwardSeedFix)
      std::cout << "  \033[1;36mBackward seed-fix: ON\033[0m\n";
   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "\n";
      }
   } else {
      std::cout << "  PID gate: NONE (fit all tracks; build proton gate from the PID plane afterwards)\n";
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
