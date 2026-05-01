/// @file run_digi_attpc.C
/// @brief Digitize the pion-in-TPC simulation: clusterize, pulse, PSA, PRA.
///
/// Reads ./data/attpcsim.root (from pi_TPC_sim.C) and writes
/// ./data/output_digi.root with AtRawEvent, AtEvent, AtPatternEvent.
///
/// Identical pipeline to 16C_pp/run_digi_attpc.C — same gas, same map,
/// same PSA/PRA settings.
///
/// Run: root -b -q run_digi_attpc.C

bool reduceFunc(AtRawEvent *evt);

void run_digi_attpc(int nEvents = 10000, float tCluster = 8.0)
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi.root";
   TString paramFile = "ATTPC.e20009_sim.par";

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

   // HELIOS-style sandbox: uniform 2x2 mm^2 square pads on a 200x200 mm^2 active
   // area centered on the beam axis. 100x100 = 10000 channels.
   auto mapping = std::make_shared<AtTpcSquareMap>(2.0, 100, 100);
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
