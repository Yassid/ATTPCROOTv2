// A/B noise-cleaning reco in ONE pass: unpack -> PSA -> {PRA on AtEventH -> AtPatternEvent (no clean)}
// AND {AtDirDeDxCleaner(AtEventH->AtEventClean) -> PRA -> AtPatternEventClean}. Same HDBSCAN config
// both branches, so the ONLY difference is the noise cleaning. Then fit both branches (fitGenfitter_ab.C).
//   root -b -q 'unpackReco_ab.C("run_0016", 2000, "/mnt/f/a1975/h5/", "/tmp/abfull/", "max")'
void unpackReco_ab(TString fileName = "run_0016", Long64_t nEvents = 2000, TString filepath = "/mnt/f/a1975/h5/",
                   TString outDir = "/tmp/abfull/", TString psaType = "max", Double_t thr = 20, int hdMcs = 20,
                   int hdMs = 8)
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer; timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString inputFile = filepath + fileName + ".h5";
   TString outputFile = outDir + fileName + "_ab_reco.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   if (gSystem->AccessPathName(inputFile.Data())) { std::cout << "ERROR: no " << inputFile << "\n"; return; }
   gSystem->Exec("mkdir -p " + outDir);
   std::cout << "IN " << inputFile << "\nOUT " << outputFile << "\n";

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

   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(false);

   AtPSA *psa = nullptr;
   if (psaType == "max") { auto p = new AtPSAMax(); p->SetThreshold(thr); psa = p; }
   else { auto p = new AtPSAMultiFit(); p->SetThreshold(thr); p->SetPeakingTime(0.720); p->SetMaxPeaks(4);
          p->SetMinSeparation(4); p->SetProminenceOnPrimary(false); psa = p; }
   psa->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a1975/D2_UKF/pad_time_correction.csv").Data());
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(false);
   psaTask->SetOutputBranch("AtEventH");

   // noise cleaning: AtEventH -> AtEventClean
   auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>(); // min_deg=1
   auto cleanTask = new AtDataCleaningTask(std::move(cleaner));
   cleanTask->SetInputBranch("AtEventH");
   cleanTask->SetOutputBranch("AtEventClean");
   cleanTask->SetPersistence(false);

   auto makeHDB = [&]() {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
      p->SetMinClusterSize(hdMcs); p->SetMinSamples(hdMs); p->SetClusterSelectionEpsilon(10.0);
      p->SetJoinMethod("mover"); p->SetMinClusterSizeJoin(15); p->SetCircleOverlapRatio(0.25);
      p->SetMotionGapTol(40); p->SetMotionAngleTol(35);
      return p;
   };
   // branch A: no clean
   auto praA = new AtPRAtask(makeHDB());
   praA->SetInputBranch("AtEventH"); praA->SetOutputBranch("AtPatternEvent"); praA->SetPersistence(true);
   // branch B: cleaned
   auto praB = new AtPRAtask(makeHDB());
   praB->SetInputBranch("AtEventClean"); praB->SetOutputBranch("AtPatternEventClean"); praB->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   run->AddTask(cleanTask);
   run->AddTask(praA);
   run->AddTask(praB);

   run->Init();
   auto numEvents = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < numEvents) numEvents = nEvents;
   std::cout << "A/B reco (unpack+" << psaType << "+2xHDBSCAN) " << numEvents << " events.\n";
   run->Run(0, numEvents);
   timer.Stop();
   std::cout << "\033[1;32mA/B reco done\033[0m -> " << outputFile << "  (" << timer.RealTime() << " s, "
             << numEvents << " evt)\n";
}
