/// @file run_ukf_only.C
/// @brief Run ONLY the UKF fitter on pre-digitized data via FairRunAna.
///
/// Reads output_digi.root (from run_digi_attpc.C) and adds the UKF fitter task.
///
/// Run: root -b -q run_ukf_only.C

void run_ukf_only(int nEvents = 10000)
{
   // Suppress verbose logging — only show errors
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi.root";
   TString outputFile = inOutDir + "output_ukf_only.root";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // UKF fitter only
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
   ukfFitter->SetZPadPlane(1000.0);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(ukfFitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(fitterTask);
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << std::endl;
   std::cout << "Macro finished succesfully." << std::endl;
   std::cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
