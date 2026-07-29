/// @file unpackFRIB_C15.C
/// @brief Unpack ONLY the FRIB DAQ group (frib/evt/*_1903) of a2091 15C HDF5 with
///        AtFRIBHDFUnpacker -> <run>_FRIB.root (AtRawEvent of the 8 auxiliary channels).
///        The ION CHAMBER is generic trace[0]; used for the 15C beam gate.
///        Fast (no PSA/PRA) — just reads the small frib group.
///
///   root -b -q 'unpackFRIB_C15.C("run_0147")'
void unpackFRIB_C15(TString fileName = "run_0147", TString outDir = "/home/yassid/a2091_C15_reco/",
                     TString filepath = "/media/yassid/Seagate Hub/ATTPC/Data/a2091/")
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = filepath + fileName + ".h5";
   TString mapDir = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = outDir + fileName + "_FRIB.root";
   TString digiParFile = dir + "/parameters/ATTPC.a2091_C15.par";
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

   // a2091 has TWO HDF5 packagings and the FRIB payload sits at a different path in each:
   //   legacy remerged : /frib/evt/evt<N>_1903                   -> AtFRIBHDFUnpacker
   //   raw merger      : /events/event_<N>/frib_physics/1903     -> AtMergerFRIBHDFUnpacker
   // Same module, same {2048, 8} shape, IC = channel 0 in both. Auto-detect, as
   // unpackReco_C15.C does for the pad data.
   std::unique_ptr<AtFRIBHDFUnpacker> unpacker;
   if (AtMergerHDFUnpacker::IsMergerFile(inputFile.Data())) {
      std::cout << "FRIB format: RAW libattpc_merger (/events/.../frib_physics/1903)"
                << " -> AtMergerFRIBHDFUnpacker" << std::endl;
      unpacker = std::make_unique<AtMergerFRIBHDFUnpacker>(fAtMapPtr);
   } else {
      std::cout << "FRIB format: legacy remerged (/frib/evt/evt<N>_1903) -> AtFRIBHDFUnpacker" << std::endl;
      unpacker = std::make_unique<AtFRIBHDFUnpacker>(fAtMapPtr);
   }
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
