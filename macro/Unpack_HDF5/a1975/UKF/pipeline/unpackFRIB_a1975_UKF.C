/// @file unpackFRIB_a1975_UKF.C
/// @brief Unpack the FRIB DAQ auxiliary stream (frib/evt group) of an a1975 run.
///
/// Reads the auxiliary HDF5 group with AtFRIBHDFUnpacker -> AtRawEvent holding
/// AtGenericTrace objects (the ion chamber, scalers, ...). Produces <run>_FRIB.root.
/// The IC amplitude per event comes from one of these generic traces; events are
/// number-ordered so they line up 1:1 with the TPC (get-stream) events.
///
/// Run: root -b -q 'unpackFRIB_a1975_UKF.C("run_0116", 5000)'

void unpackFRIB_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = -1, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = "/mnt/f/a1975/h5/" + fileName + ".h5";
   TString mapDir = dir + "/scripts/ANL2023.xml";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_FRIB.root";

   if (gSystem->AccessPathName(inputFile)) {
      printf("\033[1;31mERROR: %s not found\033[0m\n", inputFile.Data());
      return;
   }

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outputFile);
   run->SetGeomFile(dir + "/geometry/ATTPC_H1bar_geomanager.root");

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   auto unpacker = std::make_unique<AtFRIBHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(1);
   unpacker->SetBaseLineSubtraction(true);

   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(true);
   run->AddTask(unpackTask);

   run->Init();
   auto N = unpackTask->GetNumEvents();
   std::cout << "FRIB stream contains " << N << " events." << std::endl;
   if (nEvents > 0 && nEvents < N)
      N = nEvents;
   std::cout << "\033[1;32mUnpacking " << N << " FRIB-aux events.\033[0m" << std::endl;
   run->Run(0, N);

   std::cout << "\nDone. Output: " << outputFile << std::endl;
   timer.Stop();
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}
