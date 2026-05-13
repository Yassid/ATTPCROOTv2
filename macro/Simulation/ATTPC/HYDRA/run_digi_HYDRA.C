/// @file run_digi_HYDRA.C
/// @brief Digitize the HYDRA Prototype simulation: clusterize, pulse,
/// PSA, PRA.
///
/// Reads ./data/HYDRAsim<inSuffix>.root and writes
/// ./data/output_digi<outSuffix>.root. Same pipeline as pi_TPC's
/// run_digi_attpc.C with a rectangular 18×51 pad plane at 5 mm pitch
/// (AtTpcSquareMap) for the HYDRA Prototype (88 × 256 mm active).
///
/// Run: root -b -q run_digi_HYDRA.C

bool reduceFunc(AtRawEvent *evt);

void run_digi_HYDRA(int nEvents = 10000, float tCluster = 8.0,
                    bool useDriftAwareWeights = false,
                    const char *outSuffix = "",
                    bool preclusterRadiusFit = true,
                    double preclusterBin_mm = 6.0,
                    const char *inSuffix = "")
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi" + outSuffix + ".root";
   TString paramFile = "ATTPC.HYDRA_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "HYDRAsim" + inSuffix + ".root";
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

   // HYDRA Prototype MCPP: 256 × 88 mm² active area as a 128 × 44 grid
   // of 2 × 2 mm² square pads (5632 pads). X = beam direction (long),
   // Y = transverse. Lower-left vertex anchored at (x, y) = (0, 0).
   auto mapping = std::make_shared<AtTpcSquareMap>(2.0, 128, 44, 0.0, 0.0);
   mapping->GeneratePadPlane();

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetTcluster(tCluster);
   // HYDRA Prototype chamber is narrow (88 mm wide, 256 mm long). Tracks
   // span only a few cm in xy and have 10-20 hits after PSA, often in a
   // nearly-2D distribution (modest z extent — drift is 29.4 cm). The TC
   // track-finder builds 3D triplets and struggles with this flat
   // topology. Switch to the Riemann RANSAC finder, which fits circles
   // in a paraboloid transform and is more robust for nearly-planar
   // hit distributions.
   praTask->SetPRAlgorithm(3);
   praTask->SetRiemannInlierDist(15.0);  // mm — wider band to retain track tails
   praTask->SetRiemannMinHits(8);
   praTask->SetRiemannMaxTracks(1);
   praTask->SetRiemannMaxIter(500);
   praTask->SetRiemannUseSelectAndMerge(false); // single-particle events
   // Arc-walk tail recovery: after RANSAC's global inlier collection, walk
   // outward in φ from both ends of the inlier set using a sliding-window
   // local Kasa refit. Recovers track tails where the real trajectory
   // diverges from the global circle (MS, dE/dx, slight detector drift).
   praTask->SetRiemannUseArcWalkExtend(true);
   praTask->SetRiemannArcWalkWindow(30); // larger window stabilises Kasa for near-straight tracks
   praTask->SetRiemannArcWalkMaxMiss(10);
   // Replace ClusterizeSmooth3D with the arc-walk clusterizer — Smooth3D
   // produces zero clusters on these near-straight HYDRA tracks because
   // the centroid-and-step walk fails on quasi-collinear hit sets.
   praTask->SetUseArcWalk(true);
   praTask->SetTargetClusters(25);
   praTask->SetMinHitsPerCluster(3);
   praTask->SetArcWalkKNN(10);
   praTask->SetPersistence(kTRUE);
   praTask->SetUseDriftAwareWeights(useDriftAwareWeights);
   praTask->SetPreclusterRadiusFit(preclusterRadiusFit);
   praTask->SetPreclusterBin(preclusterBin_mm);
   // HYDRA: 2 mm square pads ⇒ intrinsic xy resolution = 2/√12 ≈ 0.577 mm
   praTask->SetPadResXY(2.0 / std::sqrt(12.0));

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\nMacro finished succesfully.\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}

bool reduceFunc(AtRawEvent *evt)
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
