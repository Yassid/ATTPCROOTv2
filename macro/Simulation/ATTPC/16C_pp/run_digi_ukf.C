/// @file run_digi_ukf.C
/// @brief Minimal test: adds UKF fitter to the existing digi chain.
///
/// Copy of run_digi_attpc.C with the fitter task appended.
/// Used to debug why the full pipeline hangs.

bool reduceFunc(AtRawEvent *evt);

void run_digi_ukf()
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi_ukf.root";
   TString scriptfile = "Lookup20150611.xml";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "output_digi.root"; // Pre-digitized: has AtEventH + AtPatternEvent
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();
   mapping->ParseInhibitMap("./data/inhibit.txt", AtMap::InhibitType::kTotal);

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetPersistence(kTRUE);

   // === UKF Fitter ===
   double charge = 1.602176634e-19;
   double mass = 938.272;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1));
   eloss->SetMaterial(mat);

   auto ukfFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass, std::move(eloss));
   ukfFitter->SetBField(ROOT::Math::XYZVector(0, 0, 2.85));
   ukfFitter->SetUKFParameters(1e-3, 2.0, 0.0);
   ukfFitter->SetMeasurementSigma(2.0);
   ukfFitter->SetMomentumSigmaFrac(0.3);
   ukfFitter->SetEnableEnergyStraggling(false);
   ukfFitter->SetMinClusters(10);
   // ukfFitter->SetZPadPlane(1000.0); // Disabled for debugging

   AtFitterTask *fitterTask = new AtFitterTask(std::move(ukfFitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetPersistence(kTRUE);

   // fRun->AddTask(clusterizer);
   // fRun->AddTask(pulse);
   // fRun->AddTask(psaTask);
   // fRun->AddTask(praTask);
   fRun->AddTask(fitterTask);

   fRun->Init();

   timer.Start();
   fRun->Run(0, 5);
   timer.Stop();

   std::cout << std::endl;
   std::cout << "Macro finished. Real time " << timer.RealTime() << " s" << std::endl;
}

bool reduceFunc(AtRawEvent *evt)
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
