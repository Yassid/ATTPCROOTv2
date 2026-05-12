/// @brief UKF on output_digi_perf.root, writes output_ukf_perf.root.
/// Mirrors run_ukf_only.C settings (proton, 2.85 T, H@300 torr).
/// Cluster-tangent seed enabled — not strictly needed for forward
/// stopping protons but matches the post-fix pi_TPC convention so the
/// performance plots are apples-to-apples.

void run_ukf_perf(int nEvents = 500, const char *suffix = "")
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi_perf" + suffix + ".root";
   TString outputFile = inOutDir + "output_ukf_perf" + suffix + ".root";
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

   const double e_C = 1.602176634e-19;
   double charge = e_C;
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
   ukfFitter->SetMomentumSigmaFrac(0.5);
   ukfFitter->SetEnableEnergyStraggling(false);
   ukfFitter->SetMinClusters(10);
   ukfFitter->SetNIterations(1);
   ukfFitter->SetZPadPlane(1000.0);
   // Legacy chord seed (no tangent override) — matches original
   // run_ukf_only.C. For forward stopping protons the highest-Z_digi seed
   // IS the vertex, so the chord direction is correct.
   // ukfFitter->SetUseClusterDirSeed(true);

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
   std::cout << "\nMacro finished. Real " << timer.RealTime() << " s\n";
}
