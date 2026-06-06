#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

// Unpack + PSA macro for the a1975 16C+p (proton-target) data, new-UKF workspace.
// Stage 2 of the incremental pipeline: raw HDF5 -> AtRawEvent -> AtPSAMax -> AtEventH
// (3D AtHit objects). Still no space-charge / PRA / fit. Raw events are not
// persisted here (they are large); only the reconstructed hits (AtEventH) are kept.
void unpackPSA_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = 1000, Bool_t persistRaw = false)
{
   // Load the library for unpacking and reconstruction
   gSystem->Load("libAtReconstruction.so");

   TStopwatch timer;
   timer.Start();

   TString parameterFile = "ATTPC.a1954.par";
   TString filepath = "/mnt/f/a1975/h5/"; // 16C+p raw HDF5 on the F: drive
   TString fileExt = ".h5";
   TString inputFile = filepath + fileName + fileExt;
   TString scriptfile = "ANL2023.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   TString dataDir = dir + "/macro/data/";
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());
   // _disp.root keeps the raw events too (self-contained for the event display);
   // _psa.root is the lean hits-only file for the pipeline.
   TString outputFile = persistRaw ? (fileName + "_disp.root") : (fileName + "_psa.root");
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
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   // --- Unpack ---
   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);

   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   // raw events are large; keep them only when we want a self-contained display file
   unpackTask->SetPersistence(persistRaw);

   // --- PSA (pulse-shape analysis): pads -> 3D hits ---
   auto threshold = 20;
   auto psa = new AtPSAMax();
   psa->SetThreshold(threshold);

   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(true);
   psaTask->SetOutputBranch("AtEventH");

   run->AddTask(unpackTask);
   run->AddTask(psaTask);

   std::cout << "***** Starting Init ******" << std::endl;
   run->Init();
   std::cout << "***** Ending Init ******" << std::endl;

   auto numEvents = unpackTask->GetNumEvents();
   std::cout << "Run contains " << numEvents << " events." << std::endl;
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << cGREEN << "Unpacking + PSA on " << numEvents << " events." << cNORMAL << std::endl;

   run->Run(0, numEvents);

   std::cout << std::endl << std::endl;
   std::cout << "Done unpack + PSA" << std::endl << std::endl;
   std::cout << "- Output file : " << outputFile << std::endl << std::endl;
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl << endl;
}
