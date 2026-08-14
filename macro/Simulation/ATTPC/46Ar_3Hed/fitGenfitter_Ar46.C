/// @file fitGenfitter_Ar46.C
/// @brief Fit the 46Ar(3He,d)47K simulation with EventFit::AtGenfitter, DEUTERON hypothesis.
///
/// Sibling of macro/Unpack_HDF5/a1975/D2_UKF/fitGenfitter_a1975_deuteron.C, which already fits a
/// deuteron (pdg 1000010020 -- note that file's header says "proton hypothesis", which is stale;
/// its code is a deuteron). Everything structural is copied from it. What changes is the target
/// gas, the geometry, and the field sign.
///
/// Reads <tag>_reco.root (AtPatternEvent) -- the files the accumulation already wrote. Nothing is
/// re-digitised or re-reconstructed. Writes <tag>_genfitter_d<suffix>.root with AtTrackingEvent +
/// AtPIDEvent.
///
/// B FIELD IS POSITIVE HERE, and that is not a detail. The simulation reverses the drift-z
/// handedness in digitisation, so simulated data needs bField = +2.85 where the a1975 EXPERIMENT
/// used -2.85. Measured on this very sample: reconstructed vertex z against truth gives r =
/// -1.000 with z_true + z_reco = 101.2 cm (kinematics_3Hed.C). Get the sign wrong and the fit
/// converges on mirrored tracks.
///
/// BACKWARD SEED FIX IS ON. Our deuterons run theta_lab 59 to 132 deg, so more than half are past
/// 90. On a1975 D2, without the fix ~52 % of backward tracks were seeded toward +z_lab and genfit
/// reflected them into the forward hemisphere with wrong theta and KE. The theta window is
/// therefore 50-140, not the framework default which clips at 170 but seeds forward.
///
/// NO PID GATE BY DEFAULT, ON PURPOSE. AtGenfitter owns its own AtSpyralPID (AtGenfitter.h:252)
/// and computes the gate observables itself at the class default fMinPoints = 30 -- there is no
/// SetMinPoints passthrough. A gate drawn on the mp15 plane would therefore be applied to mp30
/// observables inside the fitter, and the theta_lab 83-96 deg band that mp15 exists to recover
/// would be dropped here anyway. Since the PID is fit-independent (the same observables come out
/// of all three paths, so a gate can be re-cut without refitting), the honest order is: fit
/// everything, gate afterwards. Pass pidGate only if you have deliberately decided to live with
/// the mp30 population.
///
/// MATERIAL EFFECTS ARE OFF BY DEFAULT, inherited from a1975, and that inheritance is weaker here.
/// This gas is 8.3128e-5 g/cm3, 2.5x the H2 those macros were tuned on, and the beam loses 96 MeV
/// crossing the chamber. For an acceptance or a track-quality study kFALSE is fine and much
/// faster. For a Q-value against the proposal's 350 keV, turn it on -- and then the placeholder
/// drift velocity and the untuned gain in ATTPC.46Ar_3Hed_sim.par start to matter.
///
///   root -b -q 'fitGenfitter_Ar46.C("gs_s3001", 500)'                 // quick look
///   root -b -q 'fitGenfitter_Ar46.C("gs_s3001", -1, "/mnt/f/ar46_3hed/")'   // full sample

void fitGenfitter_Ar46(TString fileName = "gs_s3001", Long64_t nEvents = -1,
                       TString ioDir = "/mnt/f/ar46_3hed/", TString outSuffix = "", TString outDir = "",
                       Double_t bField = 2.85, Int_t minIter = 2, Int_t maxIter = 5, TString pidGate = "",
                       Double_t measSigma = 4.0, Double_t thetaMinDeg = 50.0, Double_t thetaMaxDeg = 140.0,
                       Bool_t matEffects = kFALSE, Bool_t backwardSeedFix = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   // DEUTERON ejectile for 46Ar(3He,d)47K
   const Int_t pdg = 1000010020;
   const Double_t massAmu = 2.01410177812;
   const Int_t Z = 1;

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_He3CO2_300torr_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.46Ar_3Hed_sim.par";
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter_d" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(geoManFile.Data())) {
      std::cout << "\033[1;31mERROR: " << geoManFile << " not found; run geometry/ATTPC_He3CO2_300torr.C.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfitter_Ar46 (46Ar(3He,d)47K, DEUTERON hyp) ===\033[0m\n"
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

   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << " (deuteron)\n";
   std::cout << "  gas: 3He + 5% CO2, 300 torr, 8.3128e-5 g/cm3   matEffects=" << (matEffects ? "ON" : "OFF") << "\n";
   // empty eloss file -> genfit internal; noMatEffects = !matEffects
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, /*elossFile*/ "", !matEffects,
                                                         minIter, maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);
   fitter->SetBackwardSeedFix(backwardSeedFix);
   if (backwardSeedFix)
      std::cout << "  \033[1;36mBackward seed-fix: ON\033[0m  (over half these deuterons are past 90 deg)\n";
   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "  \033[1;31m(applied at fMinPoints=30, NOT mp15)\033[0m\n";
      }
   } else {
      std::cout << "  PID gate: NONE (fit all tracks; apply the hand-drawn gate afterwards)\n";
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
