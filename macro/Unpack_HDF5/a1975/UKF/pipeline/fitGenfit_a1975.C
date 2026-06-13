/// @file fitGenfit_a1975.C
/// @brief Run ONLY the modified-GenFit fitter on a pre-reconstructed a1975 file.
///
/// Sibling of fitUKF_a1975.C: reads the SAME <fileName>_reco.root (AtPatternEvent
/// from the same PRA), but fits with AtFITTER::AtGenfit (KalmanFitterRefTrack +
/// SRIM energy loss, the modified GenFit fork) instead of the UKF. This isolates
/// the FITTER: identical PRA tracks in, genfit vs UKF kinematics out, so the Ex
/// resolution difference is purely the fitter.
///
/// Config mirrors macro/Unpack_HDF5/a1975/unpackNFit_a1975.C (the reference genfit
/// (p,d) pipeline that produced GENFIT.png): deuteron hypothesis, B=+2.85,
/// SRIM eloss proton_D2_600torr.txt, material effects OFF, merging+reclustering.
///
/// Output: <fileName>_genfit<suffix>.root with branch AtTrackingEvent holding
/// AtFittedTrack objects (GetEnergyAngles/GetExcitationEnergy/GetStats), NOT the
/// UKF kinematics format — use dres_eval_genfit.C to build the Ex spectrum.
///
/// Run: root -b -q 'fitGenfit_a1975.C("run_0106", 200, "deuteron")'
///      (geometry + par are loaded automatically from $VMCWORKDIR)

void fitGenfit_a1975(TString fileName = "run_0106", Long64_t nEvents = -1, TString particle = "deuteron",
                     Double_t magneticField = 2.85, Double_t gasMediumDensity = 0.083147, TString ioDir = "",
                     TString outSuffix = "", TString outDir = "",
                     TString elossName = "proton_D2_600torr.txt", Int_t minIter = 5, Int_t maxIter = 20,
                     Bool_t enableMerging = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";

   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfit" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m" << std::endl;
      return;
   }

   // Particle hypothesis (deuteron for (p,d); proton for (p,p'))
   Int_t pdg = 1000010020;
   Double_t mass = 2.0135532;
   Int_t Z = 1;
   if (particle == "proton") {
      pdg = 2212;
      mass = 1.00727646;
   }

   std::cout << "\033[1;33m=== fitGenfit_a1975 ===\033[0m\n";
   std::cout << "  input    : " << inputFile << "\n";
   std::cout << "  particle : " << particle << " (pdg " << pdg << ", m " << mass << ", Z " << Z << ")\n";
   std::cout << "  B        : +" << magneticField << " T   gasDensity " << gasMediumDensity << "\n";
   std::cout << "  output   : " << outputFile << "\n";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile); // needed for genfit TGeoMaterialInterface init

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   // The AtGenfit constructor inspects gGeoManager's materials, so the geometry
   // must already be loaded when the fitter is built (run->Init() loads it too
   // late). Pull the "FAIRGeom" TGeoManager from the geom file now.
   if (gROOT->FindObject("FAIRGeom") == nullptr) {
      TFile *geoFile = TFile::Open(geoManFile);
      geoFile->Get("FAIRGeom"); // becomes gGeoManager, named "FAIRGeom"
   }

   // --- GenFit fitter (mirror unpackNFit_a1975.C) ---
   Bool_t noMatEffects = 1;
   std::string elossFile = (std::string)dir.Data() + "/resources/energy_loss/" + elossName.Data();
   auto fitter = std::make_unique<AtFITTER::AtGenfit>(magneticField, 0.00001, 1000.0, elossFile, gasMediumDensity, pdg,
                                                      minIter, maxIter, noMatEffects);
   fitter->SetIonName(particle.Data());
   fitter->SetMass(mass);
   fitter->SetAtomicNumber(Z);
   fitter->SetNumFitPoints(1.0);
   fitter->SetVerbosityLevel(0);
   fitter->SetSimulationConvention(0);
   fitter->SetFitDirection(0);
   // Merging/reclustering consolidate PRA candidates and CHANGE trackIDs, which
   // breaks the downstream PID(by candidate ID)<->fit(by merged ID) match and
   // drops ~55% of gated deuterons from the analysis. Toggleable to recover them.
   fitter->EnableMerging(enableMerging);
   fitter->EnableSingleVertexTrack(enableMerging);
   fitter->EnableReclustering(enableMerging, 15.0, 7.5);

   AtFitterTaskOld *fitterTask = new AtFitterTaskOld(std::move(fitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetPersistence(kTRUE);

   // Persist the Spyral PID (brho, sqrt(dEdx), ... per track) at fit time so the
   // deuteron gate reads it back instead of re-running AtSpyralPID::Estimate() on
   // every analysis pass (the dominant analysis cost). One AtSpyralResult per PRA
   // track, parallel to GetTrackCand() order (= trackID), in the AtPIDEvent branch.
   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   run->AddTask(pidTask);
   run->AddTask(fitterTask);

   TStopwatch timer;
   timer.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   timer.Stop();

   std::cout << "\n\033[1;32mDone.\033[0m Output: " << outputFile << "  (Real " << timer.RealTime() << " s)\n";
}
