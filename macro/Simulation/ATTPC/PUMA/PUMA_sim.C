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
   puma->SetVertexXY(0., 4., 0., 4.);
   puma->SetVertexZHalfRange(22.5);
   puma->SetVertexZOffset(150.); // centre of the relocated drift volume
   puma->SetTrapRadius(10.);
   puma->SetPbarMomentum(0., 0.);
   primGen->AddGenerator(puma);
   run->SetGenerator(primGen);

   // HELIOS-style: 1 reaction per event, no separate beam-event StopTrack
   // trigger. The PUMA primaries already include the full final state.
   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   run->SetStoreTraj(kFALSE);
   run->Init();

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
