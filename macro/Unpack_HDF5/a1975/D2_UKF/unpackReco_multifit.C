// Full ATTPCROOT reconstruction (unpack -> AtPSAMultiFit -> PRA/triplclust) for one entire run,
// to regenerate the clustering with the new multipeak PSA. Output -> F:.
// Space-charge correction OFF by default (not needed; it does NOT affect z). PSA uses the Spyral
// two-point z-calibration so z matches Spyral (native CalculateZGeo compresses z ~0.55x).
// Run: root -l 'unpackReco_multifit.C("run_0300")'           // whole run
//      root -l 'unpackReco_multifit.C("run_0300", 200)'      // first 200 events
void unpackReco_multifit(TString fileName = "run_0300", Long64_t nEvents = 0, Bool_t persistRaw = true,
                         TString outDir = "/mnt/f/a1975/reco_d2/", TString filepath = "/home/yassid/spyral_d2/h5/",
                         Bool_t doSC = false, Bool_t applyTimeCorr = true, TString psaType = "multifit",
                         Double_t primSigma = 0, Double_t thr = 20, TString praType = "tc", int hdMcs = 20,
                         int hdMs = 8, Double_t fitChi2 = 0)
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer; timer.Start();

   TString parameterFile = "ATTPC.a1975_deuterium.par";
   TString inputFile = filepath + fileName + ".h5";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_" + psaType + "_reco.root";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString zlutFile = dir + "/resources/corrections/a1954/zLUT.txt";
   TString radlutFile = dir + "/resources/corrections/a1954/radLUT.txt";
   TString tralutFile = dir + "/resources/corrections/a1954/traLUT.txt";

   std::cout << "Input : " << inputFile << "\nOutput: " << outputFile << std::endl;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "ERROR: input HDF5 not found: " << inputFile << std::endl;
      return;
   }

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
   unpackTask->SetPersistence(persistRaw);

   // z via CalculateZGeo (ATTPCROOT par drift velocity) -- NOT SpyralZ. Same calibration for both PSAs.
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
      if (fitChi2 > 0) {
         p->SetFitChi2Cut(fitChi2); // fit-shape gate (float-tau + local reduced-chi2) to drop noise humps
         std::cout << "PSA   : fit-shape chi2 gate at " << fitChi2 << std::endl;
      }
      // primSigma>0: GENTLE prominence on primary at that sigma (keeps diffused far-drift pulses,
      // drops only flat baseline humps). primSigma==0: prominence off (= AtPSAMax behavior).
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
      psa->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a1975/D2_UKF/pad_time_correction.csv").Data());
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(true); // always persist full PSA cloud (AtEventH) for track viewing
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
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderHDBSCAN (mcs " << hdMcs << ", ms " << hdMs << ")" << std::endl;
   } else {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderTC>();
      p->SetClusterRadius(15.0);
      p->SetClusterDistance(7.5);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderTC (triplclust)" << std::endl;
   }
   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch(doSC ? "AtEventCorrected" : "AtEventH"); // skip SC to test z
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   if (doSC)
      run->AddTask(SCTask);
   run->AddTask(praTask);

   run->Init();
   auto numEvents = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << "Reconstructing (unpack+MultiFit+SC+triplclust) " << numEvents << " events." << std::endl;
   run->Run(0, numEvents);

   std::cout << "\nDone -> " << outputFile << std::endl;
   timer.Stop();
   std::cout << "Real time " << timer.RealTime() << " s" << std::endl;
}
