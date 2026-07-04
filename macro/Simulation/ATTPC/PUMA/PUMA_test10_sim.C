/// @file PUMA_test10_sim.C
/// @brief ATTPCROOT Geant4 transport for PUMA test channel "branch 8":
///        a single K+/K+ pair emitted back-to-back, each with fixed total
///        energy (default 0.4 GeV), momentum in the xy-plane (perpendicular
///        to the z-directed solenoid field). Clean two-track validation
///        sample — no phase space, no hexaquark.
/// Uses AtPUMAGenerator with SetChannel(10).
///
/// Run: root -b -q 'PUMA_test10_sim.C(1000)'

void PUMA_test10_sim(Int_t nEvents = 1000, Double_t testEnergy = 0.777, TString mcEngine = "TGeant4")
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

   // ---- PUMA generator, branch-10 test channel --------------------------
   // Point vertex at ~75 mm above the pad plane (middle of the readable
   // drift window). Both pions start there, back-to-back, ⊥ to B → they
   // curl into two circles at (nearly) constant z, opposite sense.
   auto *puma = new AtPUMAGenerator();
   puma->SetChannel(10);
   puma->SetTestEnergy(testEnergy); // GeV per particle
   puma->SetVertexXY(0., 0., 0., 0.);
   puma->SetVertexZHalfRange(0.);
   puma->SetVertexZOffset(75.);
   puma->SetTrapRadius(10.);
   primGen->AddGenerator(puma);
   run->SetGenerator(primGen);

   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   run->SetStoreTraj(kFALSE);
   run->Init();

   // Memory control: even a clean pi pair spawns δ-rays; keep AtMCTrack lean.
   if (auto *stack = dynamic_cast<AtStack *>(gMC->GetStack())) {
      stack->StoreSecondaries(kFALSE);
      stack->SetEnergyCut(1e-3); // GeV
      std::cout << "[PUMA_test10_sim] AtStack tuned: StoreSecondaries=false, EnergyCut=1 MeV\n";
   } else {
      std::cout << "[PUMA_test10_sim] WARNING: could not access AtStack to tune cuts.\n";
   }
   if (gMC) {
      gMC->SetCut("CUTELE", 1e-3);
      gMC->SetCut("CUTGAM", 1e-3);
      gMC->SetCut("BCUTE", 1e-3);
      gMC->SetCut("DCUTE", 1e-3);
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
   std::cout << "\nPUMA test-10 sim done (E=" << testEnergy << " GeV/particle (K+)).\n  Output: " << outFile << "\n  Real "
             << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
