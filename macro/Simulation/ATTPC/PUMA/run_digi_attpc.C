/// @file run_digi_attpc.C
/// @brief Digitize the PUMA simulation: clusterize, pulse, PSA, PRA.
///
/// Reads ./data/attpcsim.root and writes ./data/output_digi.root with
/// AtRawEvent, AtEvent, AtPatternEvent.
///
/// Mirrors 16C_pp/run_digi_attpc.C: TC (TriplClust) PRA with all defaults
/// except Tcluster. Differences vs 16C_pp are the pad map (PUMA annular
/// 16x256) and the par file (ATTPC.PUMA_sim.par).
///
/// Run: root -b -q run_digi_attpc.C

bool reduceFunc(AtRawEvent *evt);

void run_digi_attpc(int nEvents = 10000, float tCluster = 8.0)
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi.root";
   TString paramFile = "ATTPC.PUMA_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";

   TString digiParFile = dir + "/parameters/" + paramFile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // PUMA pad plane (mirrors puma-tpc-simulation/Drift/PadPlane.cpp):
   // 16 equal-area concentric rings, 256 pads per ring (4096 total), R=62.9-121.1 mm.
   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   mapping->GeneratePadPlane();

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
   praTask->SetTcluster(tCluster);
   // PUMA is annular; reaction vertex is near origin (not at high Z along beam axis as in AT-TPC).
   // Disable the primary/fragment heuristic which would otherwise reject every track.
   praTask->SetTCUseSelectAndMerge(false);
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\nMacro finished succesfully.\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}

bool reduceFunc(AtRawEvent *evt)
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
