/// @file run_digi_range.C
/// @brief Bisection helper: run digi over a [first, first+count) event
///        range with a configurable subset of tasks. Used to localize the
///        hang/leak observed around event 374.
///
/// Tasks bitmask:
///   bit 0 = AtClusterizeTask
///   bit 1 = AtPulseTask
///   bit 2 = AtPSAtask
///   bit 3 = AtPRAtask
///
/// Examples:
///   root -b -q 'run_digi_range.C(374, 1, 0x1)'   // clusterize-only on ev 374
///   root -b -q 'run_digi_range.C(370, 11, 0xF)'  // full chain, evs 370..380
///   root -b -q 'run_digi_range.C(374, 1, 0x3)'   // cluster+pulse only

void run_digi_range(int first = 374, int count = 1, unsigned mask = 0xF,
                    float tCluster = 8.0)
{
   int last = first + count;
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + Form("output_digi_range_%d_%d_m%x.root",
                                        first, last, mask);
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

   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   mapping->GeneratePadPlane();

   if (mask & 0x1) {
      auto *clusterizer = new AtClusterizeTask();
      clusterizer->SetPersistence(kFALSE);
      fRun->AddTask(clusterizer);
   }
   if (mask & 0x2) {
      auto *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
      pulse->SetPersistence(kTRUE);
      pulse->SetSaveMCInfo();
      fRun->AddTask(pulse);
   }
   if (mask & 0x4) {
      auto psa = std::make_unique<AtPSAMax>();
      psa->SetThreshold(0);
      auto *psaTask = new AtPSAtask(std::move(psa));
      psaTask->SetPersistence(kTRUE);
      fRun->AddTask(psaTask);
   }
   if (mask & 0x8) {
      auto *praTask = new AtPRAtask();
      praTask->SetTcluster(tCluster);
      praTask->SetTCUseSelectAndMerge(false);
      praTask->SetPersistence((mask & 0x10) ? kTRUE : kFALSE);
      fRun->AddTask(praTask);
   }

   fRun->Init();

   timer.Start();
   std::cout << "[run_digi_range] processing events [" << first << "," << last
             << ") with task mask 0x" << std::hex << mask << std::dec << std::endl;
   fRun->Run(first, last);
   timer.Stop();

   std::cout << "\n[run_digi_range] done. Real " << timer.RealTime()
             << " s, CPU " << timer.CpuTime() << " s\n";
}
