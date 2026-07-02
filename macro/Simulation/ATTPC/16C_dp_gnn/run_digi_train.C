// Training-data digitization WITH the additive detector-wide noise model.
// Same as run_digi_16C.C (keeps MC truth for labels: AtRawEvent SimMCPointMap +
// AtTpcPoint + AtEventH) but with baseline noise so noise hits appear and get
// labeled trackID=-1 by extract_labels.C.  Gain from par (150000).
//   root -l -b -q 'run_digi_train.C(3000, 5.0)'
void run_digi_train(Int_t nEvents = 0, Double_t noiseSigma = 5.0)
{
   TString dir = getenv("VMCWORKDIR");
   TString inOutDir = "./data/";
   TString mcFile = inOutDir + "attpcsim.root";
   TString outputFile = inOutDir + "output_digi_train.root";
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

   auto pulseObj = std::make_shared<AtPulse>(mapping);
   pulseObj->SetNoiseSigma(noiseSigma);
   pulseObj->SetNoiseAllPads(kTRUE);
   pulseObj->SetNoiseKeepThreshold(19.5); // drop noise pads below ~PSA thr=20 (global max >= windowed max => no hit lost)
   AtPulseTask *pulse = new AtPulseTask(pulseObj);
   pulse->SetPersistence(kTRUE);             // needed for SimMCPointMap (labels)
   pulse->SetSaveMCInfo();                   // pad -> mcPointID map
   pulse->SetPersistenceAtTpcPoint(kTRUE);   // AtMCPoint truth (mcPointID -> trackID)

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
   printf("\033[1;32mTrain digi (noise s=%.1f) done\033[0m -> %s  (Real %.1fs)\n", noiseSigma, outputFile.Data(),
          t.RealTime());
}
