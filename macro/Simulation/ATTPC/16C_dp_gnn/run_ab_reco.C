// A/B track-finding: psa AtEventH -> [optional AtDirDeDxCleaner] -> HDBSCAN PRA -> AtPatternEvent.
// Run with doClean=false and true (same HDBSCAN config) to compare track quality.
//   root -l -b -q 'run_ab_reco.C("in_psa.root","out_reco.root",true)'
void run_ab_reco(const char *psaFile, const char *outFile, bool doClean, int mcs = 20, int ms = 3)
{
   TString dir = getenv("VMCWORKDIR");
   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(psaFile));
   run->SetOutputFile(outFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open((dir + "/parameters/ATTPC.a1975_deuterium.par").Data(), "in");
   rtdb->setFirstInput(parIo1);

   TString praInput = "AtEventH";
   if (doClean) {
      auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>(); // min_deg=1
      AtDataCleaningTask *clean = new AtDataCleaningTask(std::move(cleaner));
      clean->SetInputBranch("AtEventH");
      clean->SetOutputBranch("AtEventClean");
      clean->SetPersistence(kTRUE);
      run->AddTask(clean);
      praInput = "AtEventClean";
   }

   auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
   p->SetMinClusterSize(mcs);
   p->SetMinSamples(ms);
   p->SetClusterSelectionEpsilon(10.0);
   p->SetMinClusterSizeJoin(15);
   AtPRAtask *praTask = new AtPRAtask(std::move(p));
   praTask->SetInputBranch(praInput);
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);
   run->AddTask(praTask);

   run->Init();
   TStopwatch t; t.Start();
   run->Run(0, 0);
   t.Stop();
   printf("\033[1;32mA/B reco (%s) done\033[0m -> %s  (%.1fs)\n", doClean ? "CLEAN" : "noclean", outFile, t.RealTime());
}
