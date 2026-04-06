/// @file test_arcwalk_ukf.C
/// @brief Run PRA + UKF with arc-walk or Smooth3D, then compare KE resolution.
///
/// Usage:
///   root -b -q 'test_arcwalk_ukf.C(500, false)'  // Smooth3D + UKF
///   root -b -q 'test_arcwalk_ukf.C(500, true)'   // Arc-walk + UKF

void test_arcwalk_ukf(int nEvents = 500, bool useArcWalk = true)
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString dir = getenv("VMCWORKDIR");
   TString paramFile = dir + "/parameters/ATTPC.e20009_sim.par";
   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi.root";
   TString outputFile = useArcWalk ? inOutDir + "output_arcwalk_ukf.root" : inOutDir + "output_smooth3d_ukf.root";

   std::cout << "Mode: " << (useArcWalk ? "Arc-walk" : "Smooth3D") << " + UKF" << std::endl;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(paramFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // PRA task
   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetInputBranch("AtEventH");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetTcluster(8.0);
   praTask->SetPersistence(kTRUE);

   if (useArcWalk) {
      praTask->SetUseArcWalk(true);
      praTask->SetTargetClusters(25);
      praTask->SetMinHitsPerCluster(3);
      praTask->SetArcWalkKNN(10);
   }

   // UKF fitter
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
   ukfFitter->SetMomentumSigmaFrac(0.5);
   ukfFitter->SetEnableEnergyStraggling(false);
   ukfFitter->SetMinClusters(10);
   ukfFitter->SetNIterations(1);
   ukfFitter->SetZPadPlane(1000.0);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(ukfFitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(praTask);
   fRun->AddTask(fitterTask);
   fRun->Init();

   TStopwatch timer;
   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\n" << (useArcWalk ? "Arc-walk" : "Smooth3D") << " + UKF done: "
             << timer.RealTime() << " s" << std::endl;
}
