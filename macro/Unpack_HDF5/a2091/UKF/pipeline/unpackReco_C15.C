/// @file unpackReco_C15.C
/// @brief a2091 15C(p,p') reco stage: unpack -> AtPSAMultiFit -> (clean) -> PRA.
///
/// Ported from a1975 D2_UKF/unpackReco_multifit.C for the 15C proton-target runs.
/// Produces <run>_reco.root (AtEventH point cloud + AtPatternEvent track candidates),
/// which the fast fitUKF_C15.C / GENFIT stages then fit.
///
/// Differences vs the a1975 macro (flagged):
///   * parFile     = ATTPC.a2091_C15.par  (B=2.85T, H2 300 torr, drift 1.30 cm/us)
///   * filepath    = /media/yassid/Seagate Hub/ATTPC/Data/a2091/   (external Seagate drive)
///   * applyTimeCorr = false  -- the pad_time_correction.csv is a1975-specific.
///                     No a2091 pad-time file exists yet; leave OFF until measured.
///   * doSC        = false    -- space charge does not affect z; not needed for a first pass.
///
///   root -l 'unpackReco_C15.C("run_0138")'        // whole run
///   root -l 'unpackReco_C15.C("run_0138", 500)'   // first 500 events
void unpackReco_C15(TString fileName = "run_0138", Long64_t nEvents = 0, Bool_t persistRaw = true,
                     TString outDir = "/home/yassid/a2091_C15_reco/",
                     TString filepath = "/media/yassid/Seagate Hub/ATTPC/Data/a2091/", Bool_t doSC = false,
                     Bool_t applyTimeCorr = false, TString psaType = "multifit", Double_t primSigma = 0,
                     Double_t thr = 20, TString praType = "tc", int hdMcs = 20, int hdMs = 8, Double_t fitChi2 = 0,
                     Double_t relErr = 0.1, TString parFile = "ATTPC.a2091_C15.par", Bool_t doClean = kTRUE,
                     Double_t bField = 2.85, Bool_t doPID = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   TStopwatch timer;
   timer.Start();

   TString parameterFile = parFile;
   TString inputFile = filepath + fileName + ".h5";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_reco.root";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar.root";
   TString zlutFile = dir + "/resources/corrections/a1954/zLUT.txt";
   TString radlutFile = dir + "/resources/corrections/a1954/radLUT.txt";
   TString tralutFile = dir + "/resources/corrections/a1954/traLUT.txt";

   std::cout << "Input : " << inputFile << "\nOutput: " << outputFile << std::endl;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "ERROR: input HDF5 not found: " << inputFile << std::endl;
      return;
   }
   gSystem->mkdir(outDir.Data(), kTRUE);

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   // a2091 mixes two HDF5 layouts: 28 runs are the legacy remerged format
   // (/get + /frib + /meta, read by AtHDFUnpacker) and 11 are the raw
   // libattpc_merger layout (/events/event_<N>/get_traces, read by
   // AtMergerHDFUnpacker). Pick the right reader automatically.
   std::unique_ptr<AtHDFUnpacker> unpacker;
   if (AtMergerHDFUnpacker::IsMergerFile(inputFile.Data())) {
      std::cout << "Format: RAW libattpc_merger (/events/get_traces) -> AtMergerHDFUnpacker" << std::endl;
      unpacker = std::make_unique<AtMergerHDFUnpacker>(fAtMapPtr);
   } else {
      std::cout << "Format: legacy remerged (/get/evt_data) -> AtHDFUnpacker" << std::endl;
      unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   }
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(persistRaw);

   AtPSA *psa = nullptr;
   if (psaType == "max") {
      auto p = new AtPSAMax();
      p->SetThreshold(thr);
      psa = p;
      std::cout << "PSA   : AtPSAMax" << std::endl;
   } else {
      auto p = new AtPSAMultiFit();
      p->SetThreshold(thr);
      p->SetPeakingTime(0.720);
      p->SetMaxPeaks(4);
      p->SetMinSeparation(4);
      if (fitChi2 != 0) {
         p->SetFloatPeakingTime(true);
         p->SetFitRelErr(relErr);
         if (fitChi2 > 0)
            p->SetFitChi2Cut(fitChi2);
         std::cout << "PSA   : float-tau, chi2 " << (fitChi2 > 0 ? "gate" : "diag(no gate)") << std::endl;
      }
      if (primSigma > 0) {
         p->SetProminenceOnPrimary(true);
         p->SetSeedSigma(primSigma);
         std::cout << "PSA   : prominence-on-primary at " << primSigma << " sigma" << std::endl;
      } else {
         p->SetProminenceOnPrimary(false);
      }
      psa = p;
      std::cout << "PSA   : AtPSAMultiFit" << std::endl;
   }
   if (applyTimeCorr)
      psa->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a2091/UKF/pipeline/pad_time_correction.csv").Data());
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(true);
   psaTask->SetOutputBranch("AtEventH");

   auto SCModel = std::make_unique<AtEDistortionModel>();
   SCModel->SetCorrectionMaps(zlutFile.Data(), radlutFile.Data(), tralutFile.Data());
   auto SCTask = new AtSpaceChargeCorrectionTask(std::move(SCModel));
   SCTask->SetInputBranch("AtEventH");

   std::unique_ptr<AtPATTERN::AtPRA> praAlgo;
   if (praType == "hdbscan") {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
      p->SetMinClusterSize(hdMcs);
      p->SetMinSamples(hdMs);
      p->SetClusterSelectionEpsilon(10.0);
      p->SetJoinMethod("mover");
      p->SetMinClusterSizeJoin(15);
      p->SetCircleOverlapRatio(0.25);
      p->SetMotionGapTol(40);
      p->SetMotionAngleTol(35);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderHDBSCAN mover" << std::endl;
   } else {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderTC>();
      p->SetClusterRadius(15.0);
      p->SetClusterDistance(7.5);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderTC (triplclust)" << std::endl;
   }

   std::string praSrc = doSC ? "AtEventCorrected" : "AtEventH";
   AtDataCleaningTask *cleanTask = nullptr;
   if (doClean) {
      auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>();
      cleanTask = new AtDataCleaningTask(std::move(cleaner));
      cleanTask->SetInputBranch(praSrc);
      cleanTask->SetOutputBranch("AtEventClean");
      cleanTask->SetPersistence(false);
      praSrc = "AtEventClean";
      std::cout << "CLEAN : AtDirDeDxCleaner (noise removal before PRA)" << std::endl;
   }

   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch(praSrc);
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);

   // PID in the chain, not in the macros. AtSpyralPID is expensive (per-track first-arc
   // circle fit + smoothing spline with a dense solve), and every downstream macro
   // -- pid_C15.C, gate_events_C15.C, ex_C15.C, brho_theta_C15.C, ic_gated_pid_C15.C --
   // was recomputing it from scratch, each one re-reading the multi-GB reco. Computing it
   // here costs almost nothing (the traces are already in memory) and persists an
   // AtPIDEvent branch carrying trackID/direction/dE/vertex/center, so a gate can be
   // APPLIED per track later instead of being recomputed to find the tracks again.
   AtPIDTask *pidTask = nullptr;
   if (doPID) {
      pidTask = new AtPIDTask();
      pidTask->SetBField(bField);
      pidTask->SetInputBranch("AtPatternEvent");
      pidTask->SetOutputBranch("AtPIDEvent");
      pidTask->SetPersistence(true);
      std::cout << "PID   : AtPIDTask (Spyral + classic) -> persistent AtPIDEvent branch, B=" << bField << " T"
                << std::endl;
   }

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   if (doSC)
      run->AddTask(SCTask);
   if (doClean)
      run->AddTask(cleanTask);
   run->AddTask(praTask);
   if (pidTask)
      run->AddTask(pidTask);

   run->Init();
   auto numEvents = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << "Reconstructing " << numEvents << " events." << std::endl;
   run->Run(0, numEvents);

   std::cout << "\nDone -> " << outputFile << std::endl;
   timer.Stop();
   std::cout << "Real time " << timer.RealTime() << " s" << std::endl;
}
