/// @file HYDRA_sim.C
/// @brief Generate charged pions inside the HYDRA Prototype (GLAD-TPC)
/// active volume.
///
/// Prototype geometry under ATTPCROOT's drift-along-z convention:
///   active volume 25.6 × 8.8 × 29.4 cm (x, y, z) — pad plane is (x,y),
///   drift length 29.4 cm along +z. Beam axis is HYDRA-z → ATTPCROOT-x
///   (the 256 mm long axis), transverse HYDRA-x → ATTPCROOT-y (88 mm).
///
/// Single-pion events (multiplicity 1 by default): vertex outside the
/// TPC upstream, pion enters through the entrance window with a small
/// angular spread around the long axis. Pass multiplicity > 1 for
/// hyperon-decay-like fan topologies.
///
/// Run: root -b -q 'HYDRA_sim.C(1000, +1)'   for pi+
///      root -b -q 'HYDRA_sim.C(1000, -1, 600., 750.)'   ~800 MeV/c pi-

void HYDRA_sim(Int_t nEvents = 1000, Int_t pionSign = +1, Double_t keMin_MeV = 600.0,
               Double_t keMax_MeV = 750.0, TString mcEngine = "TGeant4",
               Double_t Bz_T = 2.0, const char *outSuffix = "",
               Int_t multiplicity = 1)
{
   TString outFile = TString("./data/HYDRAsim") + outSuffix + ".root";
   TString parFile = TString("./data/HYDRApar") + outSuffix + ".root";

   TStopwatch timer;
   timer.Start();

   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);
   run->SetOutputFile(outFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   run->SetMaterials("media.geo");

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   AtTpc *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("ATTPC_HYDRA_Prototype_P10_1bar.root");
   run->AddModule(ATTPC);

   // Magnetic field along z (= HYDRA-y drift axis). Active region with
   // lower-left at (0, 0, 0) spans x ∈ [0, 25.6] (beam), y ∈ [0, 8.8]
   // (transverse), z ∈ [0, 29.4] (drift). Region encloses active + margin.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., Bz_T * 10.); // kG
   fMagField->SetFieldRegion(-2, 30, -2, 12, -2, 32); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   const Int_t pdg = (pionSign >= 0) ? 211 : -211;
   auto boxGen = new FairBoxGenerator(pdg, multiplicity);

   // HYDRA Prototype target: OUTSIDE the TPC, upstream of the entrance
   // window. Active region's lower-left corner is at (0, 0, 0); chamber
   // spans x ∈ [0, 25.6] (beam), y ∈ [0, 8.8] (transverse). Place the
   // target at x = -40 cm (~28 cm upstream of the entrance face),
   // y = active_y/2 = 4.4 cm (centered transverse), z = active_z/2 =
   // 14.7 cm (mid-drift).
   boxGen->SetXYZ(-40.0, 4.4, 14.7);

   // Track direction: pions emerge going downstream (+x) at θ ≈ 90°
   // (momentum in the pad plane) with a tight ±3° spread in φ around
   // +x. At 256 mm length, ±3° gives lateral excursion ≤ 13 mm — well
   // within the ±44 mm transverse half-width.
   boxGen->SetThetaRange(89.5, 90.5);
   boxGen->SetPhiRange(-3., 3.); // ±3° around +x
   boxGen->SetEkinRange(keMin_MeV * 1e-3, keMax_MeV * 1e-3); // GeV

   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);

   // Don't let AtTpc::ProcessHits StopTrack the primary on first energy
   // deposit (it interprets that as the beam reaching its reaction vertex).
   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   run->SetStoreTraj(kFALSE);
   run->Init();

   // Suppress secondary production so pad hits only contain the primary
   // pion's energy deposits. AtStack drops secondaries below 1 GeV
   // (none of the relevant secondaries exceed that), and gMC->SetCut
   // raises Geant4's production cuts so δ-rays / brems / gammas fold
   // into the primary's continuous energy loss instead of becoming
   // separate tracks.
   if (auto *stack = dynamic_cast<AtStack *>(gMC->GetStack())) {
      stack->StoreSecondaries(kFALSE);
      stack->SetEnergyCut(1.0); // GeV
      std::cout << "[HYDRA_sim] AtStack: StoreSecondaries=false, EnergyCut=1 GeV\n";
   }
   if (gMC) {
      gMC->SetCut("CUTELE", 1.0); // electrons
      gMC->SetCut("CUTGAM", 1.0); // gammas
      gMC->SetCut("BCUTE", 1.0);
      gMC->SetCut("BCUTM", 1.0);
      gMC->SetCut("DCUTE", 1.0);
      gMC->SetCut("DCUTM", 1.0);
      gMC->SetCut("PPCUTM", 1.0);
   }

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);
   run->CreateGeometryFile("./data/geofile_HYDRA_full.root");

   timer.Stop();
   std::cout << "\nMacro finished succesfully.\n";
   std::cout << "Output file:    " << outFile << "\n";
   std::cout << "Parameter file: " << parFile << "\n";
   std::cout << "PDG:            " << pdg << "  (KE: " << keMin_MeV << "-" << keMax_MeV << " MeV)\n";
   std::cout << "B field:        " << Bz_T << " T along z\n";
   std::cout << "Vertex:         (-40, 4.4, 14.7) cm — ~28 cm upstream of TPC face\n";
   std::cout << "Multiplicity:   " << multiplicity << " pions per event\n";
   std::cout << "Momentum in:    (x,y) pad plane, φ ∈ [-3°, 3°] (±3° around +x)\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n\n";
}
