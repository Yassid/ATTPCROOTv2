/// @file reco_pra_algo.C
/// @brief Run the PUMA reco chain (clusterize -> pulse -> PSA[max]) then a SELECTABLE
///        pattern-recognition algorithm, for an apples-to-apples PRA comparison
///        (same PSA hits). algo = "smooth3d" (default clustering), "arcwalk"
///        (gap-immune arc-walk clustering), or "riemann" (Riemann circle fit).
///        Writes AtPatternEvent to the output file for pra_efficiency.C.
/// Run: root -b -q 'reco_pra_algo.C("smooth3d","data/attpcsim.root","data/pra_smooth3d.root")'
void reco_pra_algo(TString algo = "smooth3d", TString mcFile = "data/attpcsim.root",
                   TString out = "data/pra_algo.root", float tCluster = 8.0)
{
   TString dir = getenv("VMCWORKDIR");
   auto *run = new FairRunAna();
   auto *source = new FairFileSource(mcFile);
   run->SetSource(source);
   run->SetOutputFile(out);
   run->SetSink(new FairRootFileSink(out));

   auto *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/ATTPC.PUMA_sim.par").Data(), "in");
   rtdb->setFirstInput(parIo);

   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   mapping->GeneratePadPlane();

   auto *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);
   auto atPulse = std::make_shared<AtPulse>(mapping);
   auto *pulse = new AtPulseTask(atPulse);
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();
   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);
   auto *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   // ---- selectable pattern recognition ----
   FairTask *praTask = nullptr;
   if (algo == "ransac") {
      // standalone sample-consensus (RANSAC) with a 2D-circle model for the curved tracks
      auto method = std::make_unique<SampleConsensus::AtSampleConsensus>(
         SampleConsensus::Estimators::kRANSAC, AtPatterns::PatternType::kCircle2D, RandomSample::SampleMethod::kUniform);
      method->SetDistanceThreshold(6.0);
      method->SetMinHitsPattern(8);
      method->SetChargeThreshold(-1);
      method->SetNumIterations(500);
      auto *sac = new AtSampleConsensusTask(std::move(method));
      sac->SetInputBranch("AtEventH");
      sac->SetPersistence(kTRUE);
      praTask = sac;
   } else {
      auto *pra = new AtPRAtask();
      pra->SetTcluster(tCluster);
      pra->SetTCUseSelectAndMerge(false); // annular geometry: primary/fragment heuristic off
      pra->SetChargeFromCenter(true);     // robust charge sign (angular sweep about the circle centre)
      pra->SetMinNumHits(6);
      if (algo == "arcwalk") { pra->SetUseArcWalk(true); pra->SetTargetClusters(8); }
      else if (algo == "riemann") { pra->SetPRAlgorithm(3); }
      else { pra->SetUseArcWalk(false); } // smooth3d
      pra->SetPersistence(kTRUE);
      praTask = pra;
   }

   run->AddTask(clusterizer);
   run->AddTask(pulse);
   run->AddTask(psaTask);
   run->AddTask(praTask);
   run->Init();
   run->Run(0, 0);
   printf("PRA_ALGO_DONE algo=%s out=%s\n", algo.Data(), out.Data());
}
