// Event display for the a1975 16C+p new-UKF workspace.
// Reads a self-contained file produced by unpackPSA_a1975_UKF(..., persistRaw=true)
// that contains BOTH AtRawEvent (pad traces) and AtEventH (3D hits, already
// reconstructed with the correct drift velocity). No on-the-fly PSA is run, so no
// parameter file is needed. Use the sidebar "Branch" control to pick which event /
// raw branch to draw.
//
// Usage:  root -l 'run_eve_UKF.C("run_0116_disp")'
void run_eve_UKF(TString InputDataFileName = "run_0116_disp")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtEventDisplay.so");

   TString InputDataFile = InputDataFileName + ".root";
   TString OutputDataFile = InputDataFileName + ".reco_display.root";
   std::cout << "Opening: " << InputDataFile << std::endl;

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = "ATTPC_H1bar_geomanager.root"; // matches the a1975 16C+p H 1 bar gas
   TString mapFile = "ANL2023.xml";

   TString GeoDataPath = dir + "/geometry/" + geoFile;
   TString mapDir = dir + "/scripts/" + mapFile;

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
   FairFileSource *source = new FairFileSource(InputDataFile);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   rtdb->setFirstInput(parIo1);

   auto fMap = std::make_shared<AtTpcMap>();
   fMap->ParseXMLMap(mapDir.Data());
   auto eveMan = new AtViewerManager(fMap);

   // 3D hit view (draws the AtEventH branch already in the file)
   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100); // max number of multihits shown

   // Pad-trace view: raw ADC + baseline-subtracted ADC + reconstructed hits
   auto tabPad = std::make_unique<AtTabPad>(2, 1);
   tabPad->DrawRawADC(0, 0);
   tabPad->DrawADC(0, 1);

   eveMan->AddTab(std::move(tabMain));
   eveMan->AddTab(std::move(tabPad));

   eveMan->Init();

   std::cout << "Finished init" << std::endl;
}
