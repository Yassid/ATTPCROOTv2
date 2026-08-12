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
/// @param useSC  apply the space-charge correction (default, and what the production used). With
///               kFALSE the SC task is not added and PRA reads AtEventH directly. This exists for
///               one diagnostic: the SC correction shifts every hit's z by a position-dependent
///               amount, which turns the time-bucket-quantised z into a continuous one. AtSpyralPID
///               collapses z-ties with EXACT equality, so on corrected data nothing merges, knot
///               gaps of ~1e-9 mm survive, and AtSmoothingSpline::Fit hits a singular pivot -- 736
///               of 10480 tracks in runs 0106-0108 (fail code 3), and they are the longest and
///               highest-charge tracks in the sample. Running with useSC = kFALSE is the control.
void unpackReco_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = 1000, Bool_t persistRaw = true,
                          TString outDir = "", Bool_t useSC = kTRUE)
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
   // PRA via dependency injection: AtPRAtask takes the algorithm as a unique_ptr
   // (mirrors AtFitterTask). AtTrackFinderTC already carries the standard TC defaults
   // the legacy AtPRAtask pushed, so only the cluster radius/distance are set here.
   // Behaviour-identical to the legacy path (validated by testPRAinject_a1975.C).
   auto praAlgo = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   praAlgo->SetClusterRadius(clusterRadius);
   praAlgo->SetClusterDistance(clusterDistance);
   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch(useSC ? "AtEventCorrected" : "AtEventH");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   if (useSC)
      run->AddTask(SCTask);
   else
      std::cout << "*** SPACE-CHARGE CORRECTION OFF -- PRA reads AtEventH (diagnostic) ***" << std::endl;
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
