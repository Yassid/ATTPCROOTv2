/// @file unpackFRIB_C14.C
/// @brief Unpack ONLY the FRIB DAQ group (frib/evt/*_1903) of a1954 14C HDF5 with
///        AtFRIBHDFUnpacker -> <run>_FRIB.root (AtRawEvent of the 8 auxiliary channels).
///        The ION CHAMBER is generic trace[0]; used for the 14C beam gate.
///        Fast (no PSA/PRA) — just reads the small frib group.
///
///   root -b -q 'unpackFRIB_C14.C("run_0147")'
void unpackFRIB_C14(TString fileName = "run_0147", TString outDir = "/home/yassid/a1954_C14_reco/",
                     TString filepath = "/mnt/f/a1954_remerged/")
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = filepath + fileName + ".h5";
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_FRIB.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954_C14.par";
   TString geoManFile = dir + "/geometry/ATTPC_H1bar.root";

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

   auto unpacker = std::make_unique<AtFRIBHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(1);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->Init();
   auto numEvents = unpackTask->GetNumEvents();
   std::cout << "Unpacking FRIB for " << numEvents << " events -> " << outputFile << std::endl;
   run->Run(0, numEvents);

   timer.Stop();
   std::cout << "\nDone -> " << outputFile << "  (Real " << timer.RealTime() << " s)" << std::endl;
}
