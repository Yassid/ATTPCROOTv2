/// @file test_arcwalk.C
/// @brief Run PRA with arc-walk or Smooth3D clustering.
///
/// Usage:
///   root -b -q 'test_arcwalk.C(500, false)'   // Smooth3D
///   root -b -q 'test_arcwalk.C(500, true)'     // Arc-walk
///
/// Then compare with: root -b -q test_arcwalk_compare.C

void test_arcwalk(int nEvents = 500, bool useArcWalk = true)
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString dir = getenv("VMCWORKDIR");
   TString paramFile = dir + "/parameters/ATTPC.e20009_sim.par";
   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi.root";
   TString outputFile = useArcWalk ? inOutDir + "output_arcwalk.root" : inOutDir + "output_smooth3d.root";

   std::cout << "Mode: " << (useArcWalk ? "Arc-walk" : "Smooth3D") << std::endl;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(paramFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

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

   fRun->AddTask(praTask);
   fRun->Init();

   TStopwatch timer;
   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\n" << (useArcWalk ? "Arc-walk" : "Smooth3D") << " done: "
             << timer.RealTime() << " s real, " << timer.CpuTime() << " s CPU" << std::endl;
}
