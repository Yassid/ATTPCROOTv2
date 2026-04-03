/// @file run_digi_noPR.C
/// @brief Digitization + PSA only (no pattern recognition).
/// Output: output_psa.root with AtRawEvent + AtEventH.
/// Then run run_pr_only.C for the pattern recognition step.

bool reduceFunc(AtRawEvent *evt);

void run_digi_noPR(int nEvents = 10000)
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_psa.root";
   TString scriptfile = "Lookup20150611.xml";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(mcFile));
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

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);

   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "Digi+PSA done: " << timer.RealTime() << " s for " << nEvents << " events" << std::endl;
}

bool reduceFunc(AtRawEvent *evt)
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
