// TEve viewer for the Dec 2014 alpha data, using the OpenKF-Claude AtViewerManager.
//
//   root -l 'eve_Dec2014_alphas.C("alpha_run_0080.root")'
//
// Tab 1 (Main)  : 3D TEve view of hits + the AT-TPC geometry, and the pad plane.
// Tab 2 (Pad)   : ADC trace of the selected pad, and its auxiliary "Q" array.
// Click a pad in the pad-plane view to load its trace into tab 2.

// withPadTab: only works on files unpacked with saveRaw=kTRUE (the ADC trace comes from
// AtRawEvent). Set it false for the bulk hits-only files in ~/dec2014_alphas_reco/.
void eve_Dec2014_alphas(TString InputDataFile = "alpha_run_0080.root",
                        TString OutputDataFile = "alpha_run_0080.display.root",
                        TString mapFile = "Lookup20141208.xml", Bool_t withPadTab = kTRUE)
{
   FairLogger *fLogger = FairLogger::GetLogger();
   fLogger->SetLogToScreen(kTRUE);
   fLogger->SetLogVerbosityLevel("MEDIUM");

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = "ATTPC_v1.1_geomanager.root";
   TString GeoDataPath = dir + "/geometry/" + geoFile;
   TString mapDir = dir + "/scripts/" + mapFile;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(InputDataFile));
   fRun->SetSink(new FairRootFileSink(OutputDataFile));
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   rtdb->setFirstInput(parIo1);

   auto fMap = std::make_shared<AtTpcMap>();
   fMap->ParseXMLMap(mapDir.Data()); // this already builds the pad plane; do not call
                                     // GeneratePadPlane() as well or it logs an error

   AtViewerManager *eveMan = new AtViewerManager(fMap);

   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100);
   eveMan->AddTab(std::move(tabMain));

   if (withPadTab) {
      auto tabPad = std::make_unique<AtTabPad>(1, 2);
      tabPad->DrawADC(0, 0);
      tabPad->DrawArrayAug("Q", 0, 1);
      eveMan->AddTab(std::move(tabPad));
   }

   eveMan->Init();
}
