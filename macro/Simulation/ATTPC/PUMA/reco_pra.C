/// @file reco_pra.C
/// @brief PUMA reconstruction through pattern recognition (no track fitting):
///        clusterize -> pulse (optional DLC dispersion) -> PSA -> PRA.
///        Reports mean pad multiplicity, PSA hits, and PRA tracks per event.
///        Fitter-free so it runs in the interpreter (the UKF fitter class is not
///        ROOT-dictionary-registered on this build).
/// Run: root -b -q 'reco_pra.C(0.0,   "base")'    // no DLC dispersion
///      root -b -q 'reco_pra.C(1.35e6,"dlc")'     // DLC @ 1.35 MOhm/sq
void reco_pra(double rDLC = 0.0, TString tag = "base", float tCluster = 8.0, double cDLC = 6.0e-13)
{
   TString dir = getenv("VMCWORKDIR");
   TString out = "./data/reco_" + tag + ".root";

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
   if (rDLC > 0)
      atPulse->SetChargeDispersionFromDLC(rDLC, cDLC, 0.5);
   auto *pulse = new AtPulseTask(atPulse);
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);
   auto *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   auto *praTask = new AtPRAtask();
   praTask->SetTcluster(tCluster);
   praTask->SetTCUseSelectAndMerge(false);
   praTask->SetPersistence(kTRUE);

   run->AddTask(clusterizer);
   run->AddTask(pulse);
   run->AddTask(psaTask);
   run->AddTask(praTask);
   run->Init();
   run->Run(0, 0);

   // Summarize reconstruction yield.
   TFile f(out);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *raw = nullptr, *evt = nullptr, *pat = nullptr;
   if (t->GetListOfBranches()->FindObject("AtRawEvent"))
      t->SetBranchAddress("AtRawEvent", &raw);
   if (t->GetListOfBranches()->FindObject("AtEventH"))
      t->SetBranchAddress("AtEventH", &evt);
   if (t->GetListOfBranches()->FindObject("AtPatternEvent"))
      t->SetBranchAddress("AtPatternEvent", &pat);
   long long n = t->GetEntries(), pads = 0, hits = 0, trks = 0, withTrk = 0;
   for (long long i = 0; i < n; i++) {
      t->GetEntry(i);
      if (raw && raw->At(0))
         pads += ((AtRawEvent *)raw->At(0))->GetNumPads();
      if (evt && evt->At(0))
         hits += ((AtEvent *)evt->At(0))->GetNumHits();
      if (pat && pat->At(0)) {
         int nt = ((AtPatternEvent *)pat->At(0))->GetTrackCand().size();
         trks += nt;
         if (nt > 0)
            withTrk++;
      }
   }
   printf(">>> RECO tag=%s  events=%lld  pads/evt=%.1f  hits/evt=%.1f  tracks/evt=%.2f  eff(>=1trk)=%.1f%%\n", tag.Data(),
          n, n ? (double)pads / n : 0, n ? (double)hits / n : 0, n ? (double)trks / n : 0,
          n ? 100.0 * withTrk / n : 0);
}
