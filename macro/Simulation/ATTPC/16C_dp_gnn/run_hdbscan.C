// Run AtTrackFinderHDBSCAN (same config as the real-data pipeline) on a point-cloud
// file (AtEventH) and write AtPatternEvent. Use on both sim (noisy) and exp to compare
// clustering behaviour under identical parameters.
//   root -l -b -q 'run_hdbscan.C("data/output_noise_s51.root","data/hdb_sim.root")'
void run_hdbscan(const char *inFile, const char *outFile, int mcs = 20, int ms = 3)
{
   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inFile));
   run->SetOutputFile(outFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
   p->SetMinClusterSize(mcs);
   p->SetMinSamples(ms);
   p->SetClusterSelectionEpsilon(10.0);
   p->SetMinClusterSizeJoin(15);
   AtPRAtask *praTask = new AtPRAtask(std::move(p));
   praTask->SetPersistence(kTRUE);

   run->AddTask(praTask);
   run->Init();
   TStopwatch t; t.Start();
   run->Run(0, 0);
   t.Stop();
   printf("\033[1;32mHDBSCAN done\033[0m %s -> %s (mcs %d ms %d, %.1fs)\n", inFile, outFile, mcs, ms, t.RealTime());
}
