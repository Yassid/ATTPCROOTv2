/// @file PUMA_sim.C
/// @brief ATTPCROOT Geant4 transport for PUMA-style events:
///        ^3He + p-bar  ->  K+ + K+ + pi- + hexaquark.
/// Uses AtPUMAGenerator (port of PUMA's GeneratePrimaries).
///
/// Run: root -b -q 'PUMA_sim.C(1000)'

void PUMA_sim(Int_t nEvents = 1000, TString mcEngine = "TGeant4")
{
   TString dir = getenv("VMCWORKDIR");
   TString outFile = "./data/attpcsim.root";
   TString parFile = "./data/attpcpar.root";

   TStopwatch timer;
   timer.Start();

   // Random seed from wall clock so re-runs give different events.
   gRandom->SetSeed(0); // 0 = use TUUID (different each run)
   std::cout << "Random seed: " << gRandom->GetSeed() << std::endl;

   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);
   run->SetOutputFile(outFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   run->SetMaterials("media.geo");

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   AtTpc *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("ATTPC_PUMA.root");
   run->AddModule(ATTPC);

   // PUMA solenoid: 4 T (= 40 kG) along z.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 40.0);
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230);
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // ---- PUMA generator -------------------------------------------------
   // Defaults match Inputs/InputAll.txt of puma-tpc-simulation:
   //   sigma_xy = 4 mm, deltaZ = 22.5 mm, p-bar at rest, He-3 target.
   // Vertex z is shifted to the centre of the relocated drift volume (z=150 mm).
   auto *puma = new AtPUMAGenerator();
   // Fixed vertex at the middle of the *readable* drift range (~75 mm above
   // the pad plane), not the geometric middle of the declared 300 mm gas.
   // At 50 MHz / NumTbs=512 / v_drift=1.5 cm/us the readout window covers
   // only ~153 mm of drift; placing the vertex at z=150 mm pushed every pulse
   // past the [20, 499] PSA search window and all hits collapsed to TB=20.
   puma->SetVertexXY(0., 0., 0., 0.);
   puma->SetVertexZHalfRange(0.);
   puma->SetVertexZOffset(75.);
   puma->SetTrapRadius(10.);
   puma->SetPbarMomentum(0., 0.);
   primGen->AddGenerator(puma);
   run->SetGenerator(primGen);

   // HELIOS-style: 1 reaction per event, no separate beam-event StopTrack
   // trigger. The PUMA primaries already include the full final state.
   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   run->SetStoreTraj(kFALSE);
   run->Init();

   // Memory control: p-bar annihilation + stopping K+ spawn thousands of
   // δ-rays per event. Without these cuts, fParticles + AtMCTrack TClonesArrays
   // grow into the millions across 1000 events and exhaust WSL memory.
   // The AtStack is created by gconfig/g4Config.C; tune it after Init().
   if (auto *stack = dynamic_cast<AtStack *>(gMC->GetStack())) {
      stack->StoreSecondaries(kFALSE); // don't persist secondaries in AtMCTrack
      stack->SetEnergyCut(1e-3);       // GeV — drop sub-MeV tracks at PushTrack
      std::cout << "[PUMA_sim] AtStack tuned: StoreSecondaries=false, EnergyCut=1 MeV\n";
   } else {
      std::cout << "[PUMA_sim] WARNING: could not access AtStack to tune cuts.\n";
   }
   // Geant4 production cuts: don't even create δ-rays/brems below 1 MeV.
   if (gMC) {
      gMC->SetCut("CUTELE", 1e-3); // electrons
      gMC->SetCut("CUTGAM", 1e-3); // gammas
      gMC->SetCut("BCUTE", 1e-3);  // electron bremsstrahlung
      gMC->SetCut("DCUTE", 1e-3);  // δ-rays from electrons
   }

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);
   run->CreateGeometryFile("./data/geofile_full.root");

   timer.Stop();
   std::cout << "\nPUMA sim done.\n  Output: " << outFile << "\n  Real " << timer.RealTime() << " s, CPU "
             << timer.CpuTime() << " s\n";
}
