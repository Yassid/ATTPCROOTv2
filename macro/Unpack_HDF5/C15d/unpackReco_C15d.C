/// @file unpackReco_C15d.C
/// @brief 15C + d  --  stage 1: raw HDF5 -> AtPatternEvent. No fitting here.
///
///   raw HDF5 -> AtRawEvent -> AtPSAMax(AtEventH) -> SC(AtEventCorrected) -> PRA(AtPatternEvent)
///
/// Self-contained: reads its own par (ATTPC.C15d_D2_300torr.par) and writes to its own
/// output area, so it shares nothing with any other workspace.
///
///   root -b -q 'unpackReco_C15d.C("run_0017", 500, true, "/home/yassid/C15d_reco/")'
///
/// The 105 usable runs are legacy-remerged HDF5 (top-level /get + /frib + /meta), which
/// AtHDFUnpacker reads directly -- no merger-format reader is needed for this set.
/// run_0090 is a TRUNCATED file and will abort; it is excluded from runs_d.txt.
///
/// persistRaw=true also keeps AtRawEvent + AtEventH so the output is self-contained for
/// the event display, at roughly 4x the file size. Use false for production.

#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

void unpackReco_C15d(TString fileName = "run_0017", Long64_t nEvents = 1000, Bool_t persistRaw = false,
                     TString outDir = "/home/yassid/C15d_reco/",
                     TString rawDir = "/media/yassid/Seagate Hub/ATTPC/Data/a1975/h5/",
                     TString parName = "ATTPC.C15d_D2_300torr.par",
                     TString geoName = "ATTPC_D300torr_v2_geomanager.root", Double_t psaThreshold = 20.0,
                     Double_t clusterRadius = 15.0, Double_t clusterDistance = 7.5,
                     // PID is computed HERE, once, and persisted gain-matched. AtSpyralPID is
                     // expensive per track (first-arc isolation, circle fit, a dense-solve smoothing
                     // spline); recomputing it in every downstream macro costs minutes per run each
                     // time and re-reads the multi-GB reco to do it.
                     // gainMatch defaults OFF: the reco persists the RAW measurement. Gain
                     // matching is a derived correction whose factors this analysis measures from
                     // its own plane (measure_gain_C15d.C), so baking it into the reco would both
                     // freeze a number that is still being determined and destroy the raw dE/dx
                     // needed to determine it. Turn it on only for a production whose table is settled.
                     Bool_t doPID = kTRUE, Bool_t gainMatch = kFALSE, TString gainTable = "",
                     // 30 = Spyral's min_total_trajectory_points for this analysis, so the two planes
                     // are cut the same way. WHATEVER VALUE A PRODUCTION PERSISTS IS THE PLANE ITS
                     // GATES MUST BE DRAWN ON -- a gate drawn at 15 does not apply to a plane cut at 30.
                     Int_t pidMinPoints = 30, Double_t pidZTieTol = 0.0, Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");

   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   if (dir.Length() == 0) {
      std::cout << cRED << "ERROR: VMCWORKDIR unset -- source build/config.sh first." << cNORMAL << std::endl;
      return;
   }
   TString scriptfile = "ANL2023.xml";
   TString mapDir = dir + "/scripts/" + scriptfile;
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());

   TString inputFile = rawDir + fileName + ".h5";
   TString outputFile = outDir + fileName + "_reco.root";
   TString digiParFile = dir + "/parameters/" + parName;
   // D2 at 300 torr. The geometry supplies genfit's MATERIAL downstream, so it has to be the
   // real target gas -- ATTPC_H1bar is ordinary hydrogen and would be wrong here even though
   // nothing complains, because with material effects off the geometry only drives navigation.
   TString geoManFile = geomDir + geoName;

   // Space-charge correction LUTs (a1954 maps, same detector and run period).
   TString zlutFile = dir + "/resources/corrections/a1954/zLUT.txt";
   TString radlutFile = dir + "/resources/corrections/a1954/radLUT.txt";
   TString tralutFile = dir + "/resources/corrections/a1954/traLUT.txt";

   std::cout << cYELLOW << "=== unpackReco_C15d (15C + d, D2 300 torr) ===" << cNORMAL << "\n"
             << "  in  : " << inputFile << "\n"
             << "  par : " << digiParFile << "\n"
             << "  geo : " << geoManFile << "\n"
             << "  out : " << outputFile << "\n";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << cRED << "ERROR: input HDF5 not found: " << inputFile << cNORMAL << std::endl;
      return;
   }
   if (gSystem->AccessPathName(digiParFile.Data())) {
      std::cout << cRED << "ERROR: par file not found: " << digiParFile << cNORMAL << std::endl;
      return;
   }
   if (gSystem->AccessPathName(geoManFile.Data())) {
      std::cout << cRED << "ERROR: geometry not found: " << geoManFile << "\n"
                << "  Generate it once with:  root -b -q " << geomDir << "ATTPC_D300torr_v2.C" << cNORMAL << std::endl;
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
   psa->SetThreshold(psaThreshold);
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(persistRaw);
   psaTask->SetOutputBranch("AtEventH");

   // --- Space-charge correction: AtEventH -> AtEventCorrected ---
   auto SCModel = std::make_unique<AtEDistortionModel>();
   SCModel->SetCorrectionMaps(zlutFile.Data(), radlutFile.Data(), tralutFile.Data());
   auto SCTask = new AtSpaceChargeCorrectionTask(std::move(SCModel));
   SCTask->SetInputBranch("AtEventH");

   // --- PRA: AtEventCorrected -> AtPatternEvent ---
   auto praAlgo = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   praAlgo->SetClusterRadius(clusterRadius);
   praAlgo->SetClusterDistance(clusterDistance);
   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch("AtEventCorrected");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   run->AddTask(SCTask);
   run->AddTask(praTask);

   // --- PID: AtPatternEvent -> AtPIDEvent, then per-run gain matching of its dE/dx ---
   if (doPID) {
      AtPIDTask *pidTask = new AtPIDTask();
      pidTask->SetInputBranch("AtPatternEvent");
      pidTask->SetOutputBranch("AtPIDEvent");
      pidTask->SetBField(bField);
      pidTask->SetMinPoints(pidMinPoints);
      pidTask->SetZTieTolerance(pidZTieTol);
      pidTask->SetPersistence(kTRUE);
      run->AddTask(pidTask);
      std::cout << "  PID  : ON (minPoints=" << pidMinPoints << ", zTieTol=" << pidZTieTol << ", B=" << bField
                << ")\n";

      if (gainMatch) {
         TString tbl = gainTable.Length() ? gainTable : (dir + "/macro/Unpack_HDF5/C15d/gainmatch_C15d.csv");
         const Int_t runNo = AtGainMatchTask::RunNumberFromName(fileName);
         if (runNo < 0) {
            std::cout << cRED << "ERROR: cannot parse a run number from '" << fileName
                      << "' -- gain matching would use the wrong factor. Aborting." << cNORMAL << std::endl;
            return;
         }
         auto *gainTask = new AtGainMatchTask(tbl.Data(), runNo);
         run->AddTask(gainTask); // AFTER AtPIDTask: it rescales that branch in place
         std::cout << "  GAIN : ON (run " << runNo << ", " << tbl << ")\n";
      } else {
         std::cout << "  GAIN : OFF (raw dE/dx persisted; gain match is applied downstream)\n";
      }
   }

   std::cout << "***** Starting Init ******" << std::endl;
   run->Init();
   std::cout << "***** Ending Init ******" << std::endl;

   auto numEvents = unpackTask->GetNumEvents();
   std::cout << "Run contains " << numEvents << " events." << std::endl;
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << cGREEN << "Reconstructing (unpack+PSA+SC+PRA) " << numEvents << " events." << cNORMAL << std::endl;

   run->Run(0, numEvents);

   timer.Stop();
   std::cout << "\n" << cGREEN << "Done." << cNORMAL << " " << outputFile << "\n"
             << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}
