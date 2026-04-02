/// @file run_reco_ukf.C
/// @brief Full reconstruction chain: digitization + PSA + pattern recognition + UKF fitting.
///
/// Runs the complete pipeline from GEANT4 MC points to UKF-fitted tracks:
///   AtClusterizeTask → AtPulseTask → AtPSAtask → AtPRAtask → AtFitterTask(UKF)
///
/// Prerequisites:
///   Run C16_pp_sim.C to generate data/attpcsim.root
///
/// Run: root -b -q run_reco_ukf.C
///      root -b -q 'run_reco_ukf.C(100)' // process 100 events

void run_reco_ukf(int nEvents = 10000)
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_reco_ukf.root";
   TString scriptfile = "Lookup20150611.xml";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   // --- Run ---
   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // --- Pad map ---
   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();
   mapping->ParseInhibitMap("./data/inhibit.txt", AtMap::InhibitType::kTotal);

   // --- Task 1: Clusterization (MC points → electron clusters) ---
   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   // --- Task 2: Pulse generation (electron clusters → pad signals) ---
   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   // --- Task 3: PSA (pad signals → hits) ---
   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   // --- Task 4: Pattern recognition (hits → tracks with clusters) ---
   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetPersistence(kTRUE);

   // --- Task 5: UKF fitter (tracks → fitted kinematics) ---
   // Proton in H2 at 300 torr
   double charge = 1.602176634e-19; // C
   double mass = 938.272;           // MeV/c^2

   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5); // g/cm^3
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1));
   eloss->SetMaterial(mat);

   auto ukfFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass, std::move(eloss));
   ukfFitter->SetBField(ROOT::Math::XYZVector(0, 0, 2.85));
   ukfFitter->SetUKFParameters(1e-3, 2.0, 0.0);
   ukfFitter->SetMeasurementSigma(2.0);
   ukfFitter->SetEnableEnergyStraggling(false); // Disable for speed in interpreted macros
   ukfFitter->SetMinClusters(5);
   ukfFitter->SetZPadPlane(1000.0); // Convert digi→lab coordinates

   AtFitterTask *fitterTask = new AtFitterTask(std::move(ukfFitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   // --- Add tasks ---
   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);
   fRun->AddTask(fitterTask);

   // --- Init and run ---
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << std::endl;
   std::cout << "Macro finished succesfully." << std::endl;
   std::cout << "Output file: " << outputFile << std::endl;

   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl;
}
