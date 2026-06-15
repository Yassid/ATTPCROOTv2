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

void fitGenfitter_a1975_deuterium(TString fileName = "run_0016", Long64_t nEvents = -1,
                                  TString ioDir = "/mnt/f/a1975/reco_d2/", TString outSuffix = "",
                                  TString outDir = "", Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5,
                                  TString pidGate = "", Double_t measSigma = 4.0, Double_t thetaMinDeg = 10.0,
                                  Double_t thetaMaxDeg = 170.0, Bool_t matEffects = kFALSE,
                                  Bool_t backwardSeedFix = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   // PROTON ejectile hypothesis for 16C(d,p)17C
   const Int_t pdg = 2212;
   const Double_t massAmu = 1.00782503207;
   const Int_t Z = 1;

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par"; // D2 gas
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter_p" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfitter_a1975_deuterium (16C(d,p)17C, PROTON hyp) ===\033[0m\n"
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

   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << " (proton)\n";
   // empty eloss file -> genfit internal; noMatEffects = !matEffects
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, /*elossFile*/ "", !matEffects,
                                                         minIter, maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg); // keep BACKWARD protons (default clips at 170)
   fitter->SetBackwardSeedFix(backwardSeedFix);      // prototype: fix backward-track seed reversal
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
