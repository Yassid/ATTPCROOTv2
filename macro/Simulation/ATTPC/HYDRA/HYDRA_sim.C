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
               Int_t multiplicity = 1, Double_t phiHalfDeg = 3.0,
               Double_t phiCenterDeg = 0.0)
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

   // Magnetic field along z (= HYDRA-y drift axis). The GLAD dipole bore
   // covers BOTH the target plane and the TPC, so the field extends well
   // beyond the chamber. The pion curves the whole way from production
   // vertex to the active region — caller must aim the pion at an off-axis
   // angle via phiCenterDeg so the curved trajectory lands inside the
   // chamber. For p = p_MeV at B = Bz_T, aim angle (HYDRA pi-, target at
   // x=-40 cm, y=4.4 cm) is φ_0 ≈ -arcsin(200/R) where R = p/(0.3·B) mm.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., Bz_T * 10.); // kG
   fMagField->SetFieldRegion(-50, 30, -15, 25, -2, 32); // cm
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
   // (momentum in the pad plane). The φ range is ±phiHalfDeg around +x;
   // pass phiHalfDeg = 0 for a delta-function-like fixed-angle beam (the
   // generator clamps to a tiny window 0 ± 0.001° to avoid degenerate
   // sampling).
   boxGen->SetThetaRange(89.5, 90.5);
   const Double_t phiHalfEff = (phiHalfDeg > 0) ? phiHalfDeg : 0.001;
   const Double_t phiLo = phiCenterDeg - phiHalfEff;
   const Double_t phiHi = phiCenterDeg + phiHalfEff;
   boxGen->SetPhiRange(phiLo, phiHi);
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
   std::cout << "Momentum in:    (x,y) pad plane, φ ∈ [" << phiLo << "°, " << phiHi << "°]\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n\n";
}
