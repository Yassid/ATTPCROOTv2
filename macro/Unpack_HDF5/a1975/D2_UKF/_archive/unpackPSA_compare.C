// Lightweight unpack -> PSA on a real a1975 event, to COMPARE point-cloud continuity
// between AtPSAMultiFit (multipeak) and AtPSAMax (single-peak). No SC/PRA (not needed here).
// Output -> F: (reco_d2). Run:
//   root -l 'unpackPSA_compare.C("run_0016", 50, "multifit")'
//   root -l 'unpackPSA_compare.C("run_0016", 50, "max")'
void unpackPSA_compare(TString fileName = "run_0016", Long64_t nEvents = 50, TString psaType = "multifit",
                       TString outDir = "/mnt/f/a1975/reco_d2/", TString filepath = "/mnt/f/a1975/h5/")
{
   gSystem->Load("libAtReconstruction.so");

   TString parameterFile = "ATTPC.a1975_deuterium.par";
   TString inputFile = filepath + fileName + ".h5";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_psa_" + psaType + ".root";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";

   std::cout << "Input : " << inputFile << "\nPSA   : " << psaType << "\nOutput: " << outputFile << std::endl;
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
   unpackTask->SetPersistence(true);

   // --- choose PSA ---
   AtPSA *psa = nullptr;
   if (psaType == "multifit") {
      auto mf = new AtPSAMultiFit();
      mf->SetThreshold(20);
      mf->SetPeakingTime(0.720);
      mf->SetMaxPeaks(4);
      mf->SetMinSeparation(4);
      mf->SetSpyralZ(560, 10, 1000); // same z frame as Spyral (for the viewer)
      mf->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a1975/D2_UKF/pad_time_correction.csv").Data());
      psa = mf;
   } else if (psaType == "mfspyral") { // AtPSAMultiFit aligned to Spyral's GET peak-selection values
      auto mf = new AtPSAMultiFit();
      mf->SetThreshold(40);       // Spyral peak_threshold
      mf->SetProminence(20);      // Spyral peak_prominence (absolute, on ALL peaks)
      mf->SetMinSeparation(50);   // Spyral peak_separation
      mf->SetPeakingTime(0.720);
      mf->SetMaxPeaks(10);        // Spyral find_peaks has no cap
      mf->SetSpyralZ(560, 10, 1000); // Spyral two-point z calib (a1975 D2): window/mm/length
      mf->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a1975/D2_UKF/pad_time_correction.csv").Data());
      psa = mf;
   } else {
      auto mx = new AtPSAMax();
      mx->SetThreshold(20);
      mx->SetSpyralZ(560, 10, 1000); // Spyral z-calibration so z is comparable to Spyral (no flip)
      mx->LoadPadTimeOffsets((dir + "/macro/Unpack_HDF5/a1975/D2_UKF/pad_time_correction.csv").Data());
      psa = mx;
   }
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(true);
   psaTask->SetOutputBranch("AtEventH");

   run->AddTask(unpackTask);
   run->AddTask(psaTask);

   run->Init();
   auto numEvents = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << "Running unpack+PSA(" << psaType << ") on " << numEvents << " events." << std::endl;
   run->Run(0, numEvents);
   std::cout << "Done -> " << outputFile << std::endl;
}
