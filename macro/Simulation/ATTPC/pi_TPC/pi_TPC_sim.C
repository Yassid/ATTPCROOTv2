/// @file pi_TPC_sim.C
/// @brief Generate charged pions inside the AT-TPC active volume.
///
/// Single-particle generator (FairBoxGenerator, PDG=+/-211). Vertex is uniform
/// inside an inscribed box of the drift tube (R=25 cm, drift_length=100 cm).
/// Direction is isotropic (theta 0-180, phi 0-360). KE drawn uniformly in
/// [keMin_MeV, keMax_MeV]. Geant4 handles dE/dx, multiple scattering and
/// (rare on TPC scales) decay-in-flight via TGeant4.
///
/// Run: root -b -q 'pi_TPC_sim.C(1000, +1)'   for pi+
///      root -b -q 'pi_TPC_sim.C(1000, -1)'   for pi-

void pi_TPC_sim(Int_t nEvents = 1000, Int_t pionSign = +1, Double_t keMin_MeV = 5.0, Double_t keMax_MeV = 50.0,
                TString mcEngine = "TGeant4", Double_t Bz_T = 2.0, const char *outSuffix = "")
{
   TString dir = getenv("VMCWORKDIR");

   TString outFile = TString("./data/attpcsim") + outSuffix + ".root";
   TString parFile = TString("./data/attpcpar") + outSuffix + ".root";

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
   ATTPC->SetGeometryFileName("ATTPC_P10_1bar.root");
   run->AddModule(ATTPC);

   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., Bz_T * 10.);               // kG (AtConstField wants kG)
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // PDG: +211 = pi+, -211 = pi-
   const Int_t pdg = (pionSign >= 0) ? 211 : -211;
   const Int_t multiplicity = 1;

   auto boxGen = new FairBoxGenerator(pdg, multiplicity);

   // Fixed vertex at the geometric centre of the drift volume
   // (R=25 cm tube, length 100 cm along Z -> centre at z=50 cm).
   boxGen->SetXYZ(0., 0., 50.);

   boxGen->SetThetaRange(5., 175.);   // isotropic (avoid pathological axial tracks)
   boxGen->SetCosTheta();             // sample uniformly in cos(theta)
   boxGen->SetPhiRange(0., 360.);

   // FairBoxGenerator::SetEkinRange takes GeV
   boxGen->SetEkinRange(keMin_MeV * 1e-3, keMax_MeV * 1e-3);

   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);

   // HELIOS-style: 1 primary particle per event, no beam/reaction. AtTpc's
   // ProcessHits otherwise StopTracks the primary as soon as fELossAcc
   // exceeds AtVertexPropagator::GetRndELoss() (default 0), interpreting
   // any energy deposit as the beam reaching its reaction vertex. Push the
   // threshold past anything achievable so the pion is tracked all the way.
   AtVertexPropagator::Instance()->SetRndELoss(1e30);

   // Visualization tracks add bulk; keep off for production.
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
   std::cout << "\nMacro finished succesfully.\n";
   std::cout << "Output file:    " << outFile << "\n";
   std::cout << "Parameter file: " << parFile << "\n";
   std::cout << "PDG:            " << pdg << "  (KE: " << keMin_MeV << "-" << keMax_MeV << " MeV)\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n\n";
}
