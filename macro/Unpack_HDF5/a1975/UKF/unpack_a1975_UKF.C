#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

// Unpack-only macro for the a1975 16C+p (proton-target) data, new-UKF workspace.
// Reads the raw HDF5 from the F: drive (NTFS, mounted at /mnt/f) and unpacks a
// limited number of events. No PSA / space-charge / PRA / fit here on purpose:
// this stage just produces the raw events so we can validate the data path and
// inspect pad traces before wiring up the UKF.
void unpack_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = 1000)
{
   // Load the library for unpacking and reconstruction
   gSystem->Load("libAtReconstruction.so");

   TStopwatch timer;
   timer.Start();

   TString parameterFile = "ATTPC.a1954.par";
   TString mappath = "";
   TString filepath = "/mnt/f/a1975/h5/"; // 16C+p raw HDF5 on the F: drive
   TString fileExt = ".h5";
   TString inputFile = filepath + fileName + fileExt;
   TString scriptfile = "ANL2023.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   TString dataDir = dir + "/macro/data/";
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());
   TString outputFile = fileName + "_unpack.root";
   TString loggerFile = dataDir + "ATTPCLog.log";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";

   std::cout << cYELLOW << "Input file : " << inputFile << cNORMAL << std::endl;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << cRED << "ERROR: input HDF5 not found: " << inputFile << cNORMAL << std::endl;
      return;
   }

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);

   // Set the parameter file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();

   std::cout << "Setting par file: " << digiParFile << std::endl;
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   std::cout << "Getting containers..." << std::endl;
   // We must get the container before initializing a run
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);

   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(true); // keep the raw events so we can inspect them

   run->AddTask(unpackTask);

   std::cout << "***** Starting Init ******" << std::endl;
   run->Init();
   std::cout << "***** Ending Init ******" << std::endl;

   // Clamp the requested number of events to what the run actually contains
   auto numEvents = unpackTask->GetNumEvents();
   std::cout << "Run contains " << numEvents << " events." << std::endl;
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << cGREEN << "Unpacking " << numEvents << " events." << cNORMAL << std::endl;

   run->Run(0, numEvents);

   std::cout << std::endl << std::endl;
   std::cout << "Done unpacking events" << std::endl << std::endl;
   std::cout << "- Output file : " << outputFile << std::endl << std::endl;
   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
   cout << endl;
   // ------------------------------------------------------------------------
}
