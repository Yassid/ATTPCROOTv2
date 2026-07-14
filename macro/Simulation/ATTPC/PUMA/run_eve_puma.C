/// @file run_eve_puma.C
/// @brief Eve-based interactive 3D event display for PUMA. Shows the annular pad
///        plane, hits + PRA tracks (3D view), the pad pulses, and the UKF & GENFIT
///        fitted trajectories.
///
///  *** NEEDS A GRAPHICAL DISPLAY (OpenGL) -- run on the desktop, NOT headless: ***
///        root 'run_eve_puma.C("data/hs_reco375.root")'
///
///  In the sidebar, pick the branches to display:
///    - event branch    -> AtEventH        (PSA hits)
///    - pattern branch   -> AtPatternEvent  (PRA track candidates)
///  The two "Fitted" tabs already point at the UKF and GENFIT tracking branches.
///  Use the entry navigator (< / >) to step through events.
void run_eve_puma(TString inputFile = "data/hs_reco375.root")
{
   TString dir = getenv("VMCWORKDIR");
   TString geoFile = dir + "/geometry/ATTPC_PUMA_geomanager.root";

   auto *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetSink(new FairRootFileSink("eve_puma_out.root"));
   fRun->SetGeomFile(geoFile);

   auto *rtdb = fRun->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/ATTPC.PUMA_sim.par").Data(), "in");
   rtdb->setFirstInput(parIo);

   // PUMA annular pad plane: 16 equal-area rings x 256 azimuthal pads (R = 62.9-121.1 mm)
   auto fMap = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   fMap->GeneratePadPlane();

   auto *eveMan = new AtViewerManager(fMap);

   // (1) 3D view: hits + PRA track candidates
   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100);
   eveMan->AddTab(std::move(tabMain));

   // (2) UKF fitted trajectories (blue), smoothed polyline, 4 T
   auto tabUKF = std::make_unique<AtTabFitted>("UKF fit", "AtTrackingEventUKF");
   tabUKF->SetDrawSmoothed(true);
   tabUKF->SetBField(4.0);
   tabUKF->SetTrackColor(kAzure + 2);
   eveMan->AddTab(std::move(tabUKF));

   // (3) GENFIT fitted trajectories (red)
   auto tabGF = std::make_unique<AtTabFitted>("GENFIT fit", "AtTrackingEventGenfit");
   tabGF->SetDrawSmoothed(true);
   tabGF->SetBField(4.0);
   tabGF->SetTrackColor(kRed + 1);
   eveMan->AddTab(std::move(tabGF));

   // (4) pad pulses: raw ADC + PSA-processed ADC
   auto tabPad = std::make_unique<AtTabPad>(2, 1, "Pulses");
   tabPad->DrawRawADC(0, 0);
   tabPad->DrawADC(0, 1);
   eveMan->AddTab(std::move(tabPad));

   eveMan->Init();
   std::cout << "PUMA event display ready. Select AtEventH / AtPatternEvent in the sidebar." << std::endl;
}
