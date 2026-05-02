/// @file run_eve.C
/// @brief Visualize the digitized PUMA events with the AtViewerManager / AtTabMain
/// event display. Reads ./data/output_digi.root and the PUMA geometry.
///
/// Run interactively (needs X / OpenGL):
///   root -l run_eve.C

void run_eve(TString InputDataFile = "./data/output_digi.root",
             TString OutputDataFile = "./data/output.reco_display.root")
{
   FairLogger *fLogger = FairLogger::GetLogger();
   fLogger->SetLogToScreen(kTRUE);
   fLogger->SetLogVerbosityLevel("MEDIUM");

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = "ATTPC_PUMA_geomanager.root";
   TString GeoDataPath = dir + "/geometry/" + geoFile;

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
   FairFileSource *source = new FairFileSource(InputDataFile);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   rtdb->setFirstInput(parIo1);

   // Same PUMA pad map used in run_digi_attpc.C — must match for the
   // event display to translate pad numbers back to (x,y) correctly.
   auto fMap = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   fMap->GeneratePadPlane();

   AtViewerManager *eveMan = new AtViewerManager(fMap);

   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100);
   eveMan->AddTab(std::move(tabMain));

   eveMan->Init();

   std::cout << "Finished init. Use the tab controls to navigate events." << std::endl;
   // eveMan->RunEvent(0);
}
