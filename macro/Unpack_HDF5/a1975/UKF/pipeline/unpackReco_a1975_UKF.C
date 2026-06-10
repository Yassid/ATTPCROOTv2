#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

// Stage 3 of the a1975 16C+p new-UKF pipeline: unpack + PSA + space-charge
// correction + PRA (pattern recognition). Produces AtPatternEvent (track
// candidates) ready for the UKF fit. Still no fit here.
//
//   raw HDF5 -> AtRawEvent -> AtPSAMax(AtEventH) -> SC(AtEventCorrected) -> PRA(AtPatternEvent)
//
// persistRaw=true also keeps AtRawEvent + AtEventH so the output is self-contained
// for the event display.
void unpackReco_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = 1000, Bool_t persistRaw = true,
                          TString outDir = "")
{
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
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());
   TString outputFile = outDir + fileName + "_reco.root";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";

   // Space-charge correction LUTs (a1954 maps, same detector/run period)
   TString zlutFile = dir + "/resources/corrections/a1954/zLUT.txt";
   TString radlutFile = dir + "/resources/corrections/a1954/radLUT.txt";
   TString tralutFile = dir + "/resources/corrections/a1954/traLUT.txt";

   std::cout << cYELLOW << "Input file : " << inputFile << cNORMAL << std::endl;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << cRED << "ERROR: input HDF5 not found: " << inputFile << cNORMAL << std::endl;
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

   // --- Unpack ---
   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(persistRaw);

   // --- PSA: pads -> 3D hits ---
   auto psa = new AtPSAMax();
   psa->SetThreshold(20);
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(persistRaw); // keep AtEventH for display when requested
   psaTask->SetOutputBranch("AtEventH");

   // --- Space-charge correction: AtEventH -> AtEventCorrected ---
   auto SCModel = std::make_unique<AtEDistortionModel>();
   SCModel->SetCorrectionMaps(zlutFile.Data(), radlutFile.Data(), tralutFile.Data());
   auto SCTask = new AtSpaceChargeCorrectionTask(std::move(SCModel));
   SCTask->SetInputBranch("AtEventH");
   // (SC task writes AtEventCorrected by default)

   // --- PRA: AtEventCorrected -> AtPatternEvent (track candidates) ---
   Double_t clusterRadius = 15.0;
   Double_t clusterDistance = 7.5;
   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetInputBranch("AtEventCorrected");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);
   praTask->SetClusterRadius(clusterRadius);
   praTask->SetClusterDistance(clusterDistance);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   run->AddTask(SCTask);
   run->AddTask(praTask);

   std::cout << "***** Starting Init ******" << std::endl;
   run->Init();
   std::cout << "***** Ending Init ******" << std::endl;

   auto numEvents = unpackTask->GetNumEvents();
   std::cout << "Run contains " << numEvents << " events." << std::endl;
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << cGREEN << "Reconstructing (unpack+PSA+SC+PRA) " << numEvents << " events." << cNORMAL << std::endl;

   run->Run(0, numEvents);

   std::cout << std::endl << "Done reconstruction" << std::endl;
   std::cout << "- Output file : " << outputFile << std::endl << std::endl;
   timer.Stop();
   cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << endl;
}
