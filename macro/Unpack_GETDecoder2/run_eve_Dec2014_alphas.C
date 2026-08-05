// TEve display for the Dec 2014 AT-TPC alpha data, on the fairroot19-port branch.
//
// This branch's AtEventDisplay is now PCL-free: the RANSAC / Hough / tracking-analysis
// overlays were removed from AtEventDrawTask (they need the PCL-only classes that
// AtReconstruction does not build here), leaving the 3D hit view, the pad plane and the
// pad-trace canvases.
//
// Reads the files produced on THIS branch by run_unpack_Dec2014_alphas.C -- the AtEvent
// schema differs between branches (ClassDef 3 here vs 6 on OpenKF-Claude), so a file
// written by one branch silently reads back with zero hits in the other.
//
//   root -l 'run_eve_Dec2014_alphas.C("/home/yassid/dec2014_alphas_reco/lowP/alpha_run_0128_hits.root")'
//
// Needs a terminal: with no tty ROOT exits as soon as the macro returns and the window goes.

// NB: the libraries are loaded by ./rootlogon.C, which ROOT runs before this macro is
// parsed. Do not add #includes or gSystem->Load calls here -- see rootlogon.C for why.
void run_eve_Dec2014_alphas(TString InputDataFile = "/home/yassid/dec2014_alphas_reco/lowP/alpha_run_0128_hits.root",
                            TString OutputDataFile = "/tmp/dec2014_alphas_display.root",
                            TString parameterFile = "ATTPC.alpha_150torr.par")
{
   FairLogger *fLogger = FairLogger::GetLogger();
   fLogger->SetLogToScreen(kTRUE);
   fLogger->SetLogVerbosityLevel("MEDIUM");

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = "ATTPC_v1.2_geomanager.root"; // the only geometry shipped on this branch
   TString GeoDataPath = dir + "/geometry/" + geoFile;

   std::cout << "Opening: " << InputDataFile << std::endl;

   // FairRoot 18.6 dropped FairRunAna::SetInputFile -- input comes from a FairSource now.
   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(InputDataFile));
   fRun->SetOutputFile(OutputDataFile);
   fRun->SetGeomFile(GeoDataPath);

   // No parameter input is set up. The FairPar*FileIo classes come out as incomplete types
   // here (only forward declarations reach cling once the Eve dictionaries are loaded), and
   // the display does not need them: AtEventDrawTask reads hit positions straight out of
   // AtEvent and builds the pad plane from AtTpcMap::GenerateAtTpc().

   AtEventManager *eveMan = new AtEventManager();

   AtEventDrawTask *eve = new AtEventDrawTask();
   eve->Set3DHitStyleBox();
   eve->SetMultiHit(100); // max multihits shown per pad
   eve->SetEventBranch("AtEventH");
   // The bulk unpack is hits-only (no AtRawEvent persisted). The raw branch is simply not
   // found in that case, which only disables the pad-trace canvases.
   eve->SetRawEventBranch("AtRawEvent");

   eveMan->AddTask(eve);
   eveMan->Init();
}
