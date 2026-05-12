/// @brief Digi pipeline writing to output_digi_perf.root (does NOT overwrite
/// the legacy 5 GB output_digi.root). Used to produce fresh PRA-fixed
/// digi data for the make_performance.C summary.

void run_digi_perf(int nEvents = 500, float tCluster = 8.0,
                   bool preclusterRadiusFit = false,
                   double preclusterBin_mm = 6.0,
                   const char *outSuffix = "")
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi_perf" + outSuffix + ".root";
   TString scriptfile = "Lookup20150611.xml";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";
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
   praTask->SetTcluster(tCluster);
   praTask->SetPersistence(kTRUE);
   praTask->SetPreclusterRadiusFit(preclusterRadiusFit);
   praTask->SetPreclusterBin(preclusterBin_mm);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();
   std::cout << "\nMacro finished. Real " << timer.RealTime() << " s\n";
}
