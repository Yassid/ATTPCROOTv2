// Preprocessing pass: direction+dE/dx noise cleaning (AtEventH -> AtEventClean).
// Insert this AtDataCleaningTask right after PSA (and after any space-charge task) in the
// reco chain; point the downstream track finder at "AtEventClean" instead of "AtEventH".
//   root -l -b -q 'run_denoise.C("in_psa.root","out.root",1)'
void run_denoise(const char *inFile = "/mnt/f/a1975/reco_d2/run_0016_psa_max.root",
                 const char *outFile = "data/denoise_test.root", int minDeg = 1)
{
   TString dir = getenv("VMCWORKDIR");
   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inFile));
   run->SetOutputFile(outFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open((dir + "/parameters/ATTPC.a1975_deuterium.par").Data(), "in");
   rtdb->setFirstInput(parIo1);

   // defaults match denoise.py: k=12, cosSeg=0.78, cosTan=0.72, rMax=32, qRatio=0.65, smoothQ=true
   auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>(12, 0.78, 0.72, 32.0, 0.65, minDeg, true);
   AtDataCleaningTask *clean = new AtDataCleaningTask(std::move(cleaner));
   clean->SetInputBranch("AtEventH");
   clean->SetOutputBranch("AtEventClean");
   clean->SetPersistence(kTRUE);

   run->AddTask(clean);
   run->Init();
   TStopwatch t; t.Start();
   run->Run(0, 0);
   t.Stop();
   printf("\033[1;32mdenoise (minDeg=%d) done\033[0m -> %s  (%.1fs)\n", minDeg, outFile, t.RealTime());
}
