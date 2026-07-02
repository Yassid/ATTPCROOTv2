// Digi with additive detector-wide baseline noise, for tuning the noise model
// against experiment. Only AtEventH (PSA hit cloud) is persisted (raw not saved)
// to keep output small. Gain comes from the par file (currently 150000).
//   root -l -b -q 'run_digi_noise.C(50, 5.0, "s5")'
void run_digi_noise(Int_t nEvents = 50, Double_t noiseSigma = 5.0, TString tag = "s5")
{
   TString dir = getenv("VMCWORKDIR");
   TString inOutDir = "./data/";
   TString mcFile = inOutDir + "attpcsim.root";
   TString outputFile = inOutDir + "output_noise_" + tag + ".root";
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
   pulseObj->SetNoiseAllPads(kTRUE); // detector-wide spray
   AtPulseTask *pulse = new AtPulseTask(pulseObj);
   pulse->SetPersistence(kFALSE); // do not save the (huge) raw event; PSA still runs in memory

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(20);
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);

   fRun->Init();
   TStopwatch t;
   t.Start();
   fRun->Run(0, nEvents);
   t.Stop();
   printf("\033[1;32mDigi(noise s=%.1f) done\033[0m -> %s  (Real %.1fs)\n", noiseSigma, outputFile.Data(), t.RealTime());
}
