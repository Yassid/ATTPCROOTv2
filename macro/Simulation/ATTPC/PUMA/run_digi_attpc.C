/// @file run_digi_attpc.C
/// @brief Digitize the PUMA simulation: clusterize, pulse, PSA, PRA, optional UKF.
///
/// Reads ./data/attpcsim.root and writes ./data/output_digi.root with
/// AtRawEvent, AtEvent, AtPatternEvent (and AtTrackingEvent + AtFitMetadata
/// when runUKF=true). Pass runUKF=true to fold the multi-hypothesis UKF fitter
/// into the same job — useful when you want a single output file the event
/// display can render without an additional run_ukf_multi step.
///
/// Run: root -b -q run_digi_attpc.C
///      root -b -q 'run_digi_attpc.C(1000, 8.0, true)'         // K+/K-/π+/π-
///      root -b -q 'run_digi_attpc.C(1000, 8.0, true, true)'   // π+/π- only

bool reduceFunc(AtRawEvent *evt);

void run_digi_attpc(int nEvents = 10000, float tCluster = 8.0,
                    bool runUKF = false, bool pionsOnly = false)
{
   TString inOutDir = "./data/";
   TString outputFile = inOutDir + "output_digi.root";
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

   // PUMA pad plane (mirrors puma-tpc-simulation/Drift/PadPlane.cpp):
   // 16 equal-area concentric rings, 256 pads per ring (4096 total), R=62.9-121.1 mm.
   auto mapping = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
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
   // PUMA is annular; reaction vertex is near origin (not at high Z along beam axis as in AT-TPC).
   // Disable the primary/fragment heuristic which would otherwise reject every track.
   praTask->SetTCUseSelectAndMerge(false);
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   // Optional inline UKF fitter — same configuration as run_ukf_multi.C.
   // Adds AtTrackingEvent + AtFitMetadata to the digi output so the event
   // display can render fitted helices without a second run.
   if (runUKF) {
      const double e_C = 1.602176634e-19;
      const double u_MeV = 931.49410372;
      const double m_K_MeV = 493.677;
      const double m_pi_MeV = 139.57039;
      // pionsOnly skips the K+/K- hypotheses — used to isolate the UKF's
      // π+/π- mass discrimination on a pure-pion sample (PiGun_sim).
      const int nHyp = pionsOnly ? 2 : 4;
      const double masses[4] = {m_K_MeV, m_K_MeV, m_pi_MeV, m_pi_MeV};
      const int signs[4] = {+1, -1, +1, -1};
      const std::string names[4] = {"K+", "K-", "pi+", "pi-"};
      const int firstHyp = pionsOnly ? 2 : 0; // skip K+/K- when pionsOnly

      auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
      for (int i = firstHyp; i < firstHyp + nHyp; ++i) {
         auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
         eloss->SetProjectile(1, 1, masses[i] / u_MeV);
         std::vector<std::tuple<int, int, int>> mat;
         mat.push_back(std::make_tuple(18, 40, 9));
         mat.push_back(std::make_tuple(6, 12, 1));
         mat.push_back(std::make_tuple(1, 1, 4));
         eloss->SetMaterial(mat);

         auto ukf = std::make_unique<EventFit::AtFitterUKF>(signs[i] * e_C, masses[i], std::move(eloss));
         ukf->SetBField(ROOT::Math::XYZVector(0, 0, 4.0));
         ukf->SetUKFParameters(1e-3, 2.0, 0.0);
         ukf->SetMeasurementSigma(1.0);
         ukf->SetMomentumSigmaFrac(0.1);
         ukf->SetEnableEnergyStraggling(true);
         ukf->SetMinClusters(2);
         ukf->SetNIterations(3);
         ukf->SetZPadPlane(0.0);
         ukf->SetBackExtrapMaxPath(250.0);
         ukf->SetUseHelixBackExtrap(true);
         ukf->SetMinClusterSpacing(1.0);
         ukf->SetTargetClusters(5);
         ukf->SetAdaptiveDistBounds(6.0, 18.0);
         ukf->SetUseClusterDirSeed(true);
         // Geometry-based, gap-immune re-clustering. Better than Smooth3D
         // for PUMA's centered-vertex topology where the Z-max seed of
         // Smooth3D lands at the wrong end of upward-going tracks.
         ukf->SetUseArcWalk(true);
         // Compensate the hardwired +8.6 mm digi-z bias in AtPSA's z-reco
         // (peakingTime · vDrift = 500 ns · 1.5 cm/μs = 7.5 mm, +1 mm
         // pulse-shape asymmetry — measured directly via digi_z_offset.C).
         ukf->SetVertexZBias(8.6);
         multi->AddHypothesis(std::move(ukf), names[i], signs[i]);
      }

      auto *fitterTask = new AtFitterTask(std::move(multi));
      fitterTask->SetInputBranch("AtPatternEvent");
      fitterTask->SetOutputBranch("AtTrackingEvent");
      fitterTask->SetFitMetadataBranch("AtFitMetadata");
      fitterTask->SetPersistence(kTRUE);
      fRun->AddTask(fitterTask);
   }

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
