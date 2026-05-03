/// @file PiGun_sim.C
/// @brief ATTPCROOT Geant4 transport for the random π+/π- gun. Same geometry
///        and B-field as PUMA_sim.C, but only pions (random multiplicity, sign,
///        direction, momentum). Use this to isolate the UKF's K-vs-π mass
///        discrimination from the PUMA topology.
///
/// Run: root -b -q 'PiGun_sim.C(1000)'

void PiGun_sim(Int_t nEvents = 1000, TString mcEngine = "TGeant4")
{
   TString dir = getenv("VMCWORKDIR");
   TString outFile = "./data/attpcsim.root";
   TString parFile = "./data/attpcpar.root";

   TStopwatch timer;
   timer.Start();

   gRandom->SetSeed(0);
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

   // ---- π+/π- gun ----------------------------------------------------------
   // Same vertex as PUMA (point source at (0, 0, 75) mm). Multiplicity 1..3
   // mirrors PUMA's ~3 charged primaries per event. Momentum 100..300 MeV/c
   // covers the PUMA pion KE range (~30..200 MeV KE) and a bit beyond so the
   // sample also exercises the higher-momentum end where mass discrimination
   // becomes harder.
   auto *gun = new AtPiGunGenerator();
   gun->SetVertexXY(0., 0., 0., 0.);
   gun->SetVertexZHalfRange(0.);
   gun->SetVertexZOffset(75.);
   gun->SetTrapRadius(10.);
   gun->SetMultiplicityRange(1, 1);  // exactly one pion per event
   gun->SetMomentumRange(0.10, 0.30);
   gun->SetCharge(0); // 0 = random ±, +1 = π+ only, -1 = π- only
   primGen->AddGenerator(gun);
   run->SetGenerator(primGen);

   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   run->SetStoreTraj(kFALSE);
   run->Init();

   // Deactivate secondaries:
   //  * AtStack drops every track pushed below 1 GeV — no MCTrack entry,
   //    no further propagation of that secondary, no MC points along its
   //    (would-be) path.
   //  * gMC->SetCut(...) raises the energy threshold below which Geant4
   //    folds secondaries into continuous loss along the primary's path
   //    instead of producing them as separate tracks.
   if (auto *stack = dynamic_cast<AtStack *>(gMC->GetStack())) {
      stack->StoreSecondaries(kFALSE);
      stack->SetEnergyCut(1.0);   // GeV
      std::cout << "[PiGun_sim] AtStack: StoreSecondaries=false, EnergyCut=1 GeV\n";
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
   run->CreateGeometryFile("./data/geofile_full.root");

   timer.Stop();
   std::cout << "\nPiGun sim done.\n  Output: " << outFile << "\n  Real " << timer.RealTime()
             << " s, CPU " << timer.CpuTime() << " s\n";
}
