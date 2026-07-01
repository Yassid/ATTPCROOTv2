// Digitize the 16C(d,p) sim into pad signals + PSA hit cloud, KEEPING MC truth track
// IDs (SaveMCInfo). Output has: AtRawEvent (with pad->mcPointID map), AtEventH (PSA
// hits), AtTpcPoint (AtMCPoint truth, mcPointID->trackID). The label extractor joins
// hit.padNum -> SimMCPointMap -> mcPointID -> AtMCPoint.trackID.
//
//   root -l -b -q 'run_digi_16C.C'
void run_digi_16C(Int_t nEvents = 0)
{
   TString dir = getenv("VMCWORKDIR");
   TString inOutDir = "./data/";
   TString mcFile = inOutDir + "attpcsim.root";
   TString outputFile = inOutDir + "output_digi.root";
   // a1975 deuterium conditions (drift velocity 1.15 cm/us) + a1975 pad map
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";
   TString mapParFile = dir + "/scripts/ANL2023.xml";

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

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();                  // pad -> mcPointID map in AtRawEvent
   pulse->SetPersistenceAtTpcPoint(kTRUE);  // keep AtMCPoint truth (mcPointID -> trackID)

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(20);
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);

   fRun->Init();
   TStopwatch t; t.Start();
   fRun->Run(0, nEvents);
   t.Stop();
   printf("\033[1;32mDigi done\033[0m -> %s  (Real %.1fs)\n", outputFile.Data(), t.RealTime());
}
