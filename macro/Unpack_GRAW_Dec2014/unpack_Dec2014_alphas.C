// Unpack the Dec 2014 AT-TPC alpha data (external drive Cris_OneD) on the OpenKF-Claude
// branch, so the modern TEve viewer (AtViewerManager) can be used on it.
//
// The .graw files are fed through symlinks named file<N>_... because AtGRAWUnpacker maps
// files to decoders by a "%i" pattern and this run has no CoBo 5 (CoBos are 0-4, 6-9).
// The real CoBo number stays in the frame header, which is what the pad mapping uses,
// so the renumbering is bookkeeping only.
//
//   root -l 'unpack_Dec2014_alphas.C("/home/yassid/dec2014_links/run_0080.txt","alpha_run_0080.root",500)'

void unpack_Dec2014_alphas(TString dataFile = "/home/yassid/dec2014_links/run_0080.txt",
                           TString outputFile = "alpha_run_0080.root", Int_t nEvents = 500,
                           Int_t numFiles = 9, TString mapFile = "Lookup20141208.xml",
                           TString parameterFile = "ATTPC.alpha.par", Bool_t saveRaw = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");

   TStopwatch timer;
   timer.Start();

   TString dir = gSystem->Getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + mapFile;
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_v1.1.root";

   FairRunAna *run = new FairRunAna();
   run->SetSink(new FairRootFileSink(outputFile));
   run->SetGeomFile(geoManFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   auto unpacker = std::make_unique<AtGRAWUnpacker>(fAtMapPtr, numFiles);
   unpacker->SetInputFileName(dataFile.Data(), "file%i_");
   unpacker->SetInitialEventID(0);
   // 2014 data carry no topology frame -- synthesise one, as the old AtDecoder2Task did
   // (mask 0xF = four AsAds per CoBo).
   unpacker->SetPseudoTopologyFrame(0xF, kFALSE);
   // Must stay true: AtPSAMax requires pedestal-subtracted pads and finds nothing without it.
   unpacker->SetSubtractFPN(true);
   unpacker->SetSaveLastCell(false);
   // NB: do NOT call SetCheckNumEvents() here. GetNumEvents() walks every file to EOF to
   // count, which leaves the decoders parked at end-of-file; GetBasicFrame can seek back but
   // GetCoboFrame cannot, so the CoBo-frame path then returns null for every event. Pass an
   // explicit event count instead (see count_events.C in the fr19port worktree).
   // Required: merge all 4 AsAd frames per CoBo into one event. Without it each event holds
   // a single AsAd (2304 pads instead of 9216) with the AsAd index cycling 0,1,2,3.
   unpacker->SetUseCoboFrame(kTRUE);

   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   // Persisting AtRawEvent costs ~1.3 MB/event (2304 pads x 512 samples). Keep it on for
   // runs you want to open in the TEve viewer's pad/ADC tab; turn it off for bulk unpacking,
   // where it is ~100x smaller and only the hits are needed.
   unpackTask->SetPersistence(saveRaw);

   auto psa = new AtPSAMax();
   psa->SetThreshold(20);

   auto psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);

   run->Init();

   if (nEvents <= 0) {
      std::cout << "Give an explicit event count (see count_events.C)." << std::endl;
      return;
   }
   std::cout << "Unpacking " << nEvents << " events." << std::endl;

   run->Run(0, nEvents);

   timer.Stop();
   std::cout << std::endl;
   std::cout << "Done. Output: " << outputFile << std::endl;
   std::cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
