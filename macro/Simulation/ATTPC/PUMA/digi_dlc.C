/// @file digi_dlc.C
/// @brief Minimal PUMA digitization (clusterize -> pulse) to compare the raw pad
///        multiplicity with and without DLC resistive charge dispersion.
/// Run: root -b -q 'digi_dlc.C(0.0,   "base")'      // no dispersion
///      root -b -q 'digi_dlc.C(1.35e6,"dlc")'       // DLC @ 1.35 MOhm/sq
void digi_dlc(double rDLC = 0.0, TString tag = "base", double cDLC = 6.0e-13)
{
   TString dir = getenv("VMCWORKDIR");
   TString out = "./data/digi_" + tag + ".root";

   auto *run = new FairRunAna();
   run->SetSource(new FairFileSource("./data/attpcsim.root"));
   run->SetOutputFile(out);

   auto *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/ATTPC.PUMA_sim.par").Data(), "in");
   rtdb->setFirstInput(parIo);

   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   mapping->GeneratePadPlane();

   auto *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   auto atPulse = std::make_shared<AtPulse>(mapping);
   if (rDLC > 0) {
      atPulse->SetChargeDispersionFromDLC(rDLC, cDLC, 0.5);
      std::cout << "[digi_dlc] DLC ON: R=" << rDLC / 1e6 << " MOhm/sq" << std::endl;
   } else {
      std::cout << "[digi_dlc] baseline (no dispersion)" << std::endl;
   }
   auto *pulse = new AtPulseTask(atPulse);
   pulse->SetPersistence(kTRUE);

   run->AddTask(clusterizer);
   run->AddTask(pulse);
   run->Init();
   run->Run(0, 0); // all events

   // Report mean pad multiplicity.
   TFile f(out);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *raw = nullptr;
   t->SetBranchAddress("AtRawEvent", &raw);
   long long tot = 0, n = t->GetEntries();
   for (long long i = 0; i < n; i++) {
      t->GetEntry(i);
      auto *ev = (AtRawEvent *)raw->At(0);
      if (ev)
         tot += ev->GetNumPads();
   }
   printf(">>> RESULT tag=%s  events=%lld  total_pads=%lld  mean_pads/event=%.1f\n", tag.Data(), n, tot,
          n ? (double)tot / n : 0.0);
}
