/// @file unpackFRIB_C15p.C
/// @brief Unpack ONLY the FRIB DAQ group (frib/evt/*_1903) of a2091 H2 (15C+p) HDF5 with
///        AtFRIBHDFUnpacker -> <run>_FRIB.root (AtRawEvent of the 8 auxiliary channels).
///        The ION CHAMBER is generic trace[0]; used for the 15C beam gate.
///        Fast (no PSA/PRA) — just reads the small frib group.
///
///   root -b -q 'unpackFRIB_C15p.C("run_0017")'
void unpackFRIB_C15p(TString fileName = "run_0017", TString outDir = "/home/yassid/a2091_C15_ic/tmp/",
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
   // This workspace's own D2 geometry, not the inherited H1bar: that file does not even exist in a
   // fresh checkout (the .root geometries are generated, not tracked), and its absence crashed this
   // macro with a segfault rather than an error. The FRIB group carries no tracks, so the geometry
   // only has to LOAD -- but it still has to exist.
   TString geoManFile = dir + "/geometry/ATTPC_H300torr_RT_geomanager.root";
   if (gSystem->AccessPathName(geoManFile.Data())) {
      std::cout << "ERROR: geometry not found: " << geoManFile
                << "\n  generate it once:  root -b -q $VMCWORKDIR/geometry/ATTPC_H300torr_RT.C\n";
      return;
   }

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
