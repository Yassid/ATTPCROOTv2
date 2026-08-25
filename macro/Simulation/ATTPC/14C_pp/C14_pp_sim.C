/// @file C14_pp_sim.C
/// @brief a1954 14C(p,p')14C simulation, matched to the experimental running conditions.
///
/// Conditions:
///     B = 2.85 T (28.5 kG),  H2 at 300 torr,  drift length 1000 mm
///     beam 14C at 161 MeV = 11.5 MeV/u
///
/// GAS PRESSURE IS 300 torr AT ROOM TEMPERATURE: rho = 3.308e-5 g/cm3, geometry
/// ATTPC_H300torr_RT.root (changed 2026-08-25; was ATTPC_H300torr.root = 3.553e-5, a 0 C value,
/// which put 7.4 % too much material in the truth and disagreed with this sim's own digi par).
/// Note that this CONTRADICTS parameters/ATTPC.a1954_C14.par, which says GasPressure 600,
/// and it contradicts the production reconstruction, which corrected energy loss with
/// 600-torr material in both fitters:
///     fitUKF_C14.C     gasDensity = 6.5e-5 g/cm3   ("H2 at 600 torr")
///     fitGenfit_C14.C  geoName    = "ATTPC_H600torr"  (rho = 6.616e-5)
/// That is ~1.9x too much material. Since the energy-loss correction scales with path
/// length, an error there is angle- and vertex-z-dependent, which is exactly the shape of
/// the two drifts seen in the data. Reconstructing this sim (truth at 300 torr) with the
/// 600-torr fitter settings is the clean test of how much of the observed g.s. offset and
/// drift that mistake accounts for.
///
/// Beam momentum: the ion generator multiplies the pz argument by A, so pz is GeV/c per
/// nucleon and the TOTAL momentum is what sets the beam energy.
///     m(14C) = 14.003242 u = 13.04394 GeV
///     KE = 161 MeV  ->  p = 2.055740 GeV/c  ->  pz = 2.055740/14 per nucleon
/// NomEnergy does NOT set the beam energy -- only pz does; AtVertexPropagator::GetBeamNomE()
/// has no consumer anywhere in the tree. But it is NOT harmless either: the same constructor
/// argument doubles as the default for maxELoss, which controls how often a reaction is
/// generated at all (see the maxELoss comment below). Pass maxELoss explicitly.
///
///   root -l 'C14_pp_sim.C(2000)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. That discrepancy
/// is repo-wide and is kept here so this sim stays comparable to its siblings -- but it
/// means the beam energy MUST be verified from MC truth rather than assumed. See
/// check_beam_C14.C, which does exactly that.

/// CM ANGULAR RANGE: default 5-120 deg, not the full 0-180. In inverse kinematics the recoil
/// proton comes out at theta_lab ~ (180 - theta_cm)/2, so theta_cm 5-120 maps to roughly
/// theta_lab 30-88 -- the region the AT-TPC actually accepts. Beyond ~120 deg CM the recoil
/// proton is too low in energy to make a reconstructable track, so generating there just
/// produces events with no track and wastes the digitization (the slow step).
void C14_pp_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 5.0, Double_t thetaMaxCM = 120.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                Double_t resEx = 0.0, UInt_t seed = 0,
                // TRANSPORT GAS. Changed 2026-08-25 from ATTPC_H300torr.root (3.553e-5) to the _RT
                // variant (3.308e-5). Both are "H2 at 300 torr"; the first is at 0 C and the second
                // at 293 K, and the a1954 gas is at ROOM TEMPERATURE. The old value put 7.4 % too
                // much material in the TRUTH, and it also disagreed with this sim's own digitisation
                // par (ATTPC.a1954_C14_sim.par: Density 0.0331 mg/cm3 = 3.31e-5) -- transport and
                // digitisation were modelling different gases. Unlike a fitter geometry with
                // material effects off, this is NEVER inert: it sets the generated energy loss.
                TString geoFile = "ATTPC_H300torr_RT.root")
{
   // RNG seed. There was NO seeding here, so parallel jobs would have produced byte-identical
   // events and any "added statistics" would have been a copy of the same sample. seed = 0 keeps
   // ROOT's default (time-based) behaviour; pass a distinct value per parallel job.
   // The same pattern as macro/Simulation/ATTPC/15C_d/C15_dp_sim.C:6.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed. gRandom->GetSeed() returns TRandom3's internal state counter
   // (always 624 right after SetSeed), not the seed -- useless for verifying a parallel run.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar.root";

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

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName(geoFile); // default ATTPC_H300torr_RT.root, rho = 3.308e-5 g/cm3
   std::cout << "\033[1;33m[C14_pp_sim] transport gas geometry: " << geoFile << "\033[0m" << std::endl;
   run->AddModule(ATTPC);

   // -----   Magnetic field : 2.85 T, as in ATTPC.a1954_C14.par   -----------
   // SIGN: the a1954 DATA corresponds to -2.85 T. Generating at +28.5 kG gave the sim the
   // opposite sense of rotation, and AtSpyralPID::Direction() infers forward/backward purely
   // from that sense (it never looks at z) while `polar` comes from d(rho)/dz -- so its
   // consistency check rejected 82 % of sim tracks (only 10.8 % valid, vs 100 % on data).
   // Mirroring a transverse coordinate, i.e. flipping the rotation sense alone, restored
   // 96.8 %. Generate in the DATA convention so nothing downstream needs a sign special-case.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 14C at 161 MeV (11.5 MeV/u)   ---------------------------
   Int_t z = 6;  // Atomic number
   Int_t a = 14; // Mass number
   Int_t q = 0;  // Charge state
   Int_t m = 1;  // Multiplicity
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = 2.055740 / a; // GeV/c per nucleon -> 161.0 MeV total KE
   Double_t BExcEner = 0.0;
   Double_t Bmass = 14.003242; // amu (repo convention -- see the note above)
   Double_t NomEnergy = 161.0; // MeV

   // maxELoss CONTROLS THE REACTION RATE -- it is not cosmetic.
   // AtTPCIonGenerator does  fMaxEnLoss = (eLoss < 0 ? ener : eLoss)  and then samples
   //     Er = gRandom->Uniform(0, fMaxEnLoss);  SetRndELoss(Er)
   // while AtTpc::reactionOccursHere() fires only once  fELossAcc > GetRndELoss().
   // A 14C beam at 161 MeV loses only ~12 MeV crossing the whole 1 m of H2 at 300 torr, so
   // leaving eLoss at its -1 default (=> threshold uniform in [0,161]) means the threshold
   // exceeds the available energy loss ~92% of the time and NO reaction is generated:
   // measured 300 reactions in 8000 events (7.5%), matching 12/161.
   // Setting it to the actual traversal loss makes essentially every reaction event fire,
   // with the vertex still uniform in z.
   Double_t maxELoss = 12.0; // MeV, ~ the 14C energy loss across the full drift length

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 14C(p,p')14C   -------------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, scattered, recoil -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 14C ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(14.003242); // uma
   ExE.push_back(BExcEner);

   // ---- Target : p ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.0078250322); // uma
   ExE.push_back(0.0);

   // ---- Scattered : 14C ----
   Zp.push_back(6);
   Ap.push_back(14);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(14.003242); // uma
   // Excitation of the RESIDUAL 14C. 0 = elastic; 6.094 = the first excited state (1-).
   // Inelastic changes the two-body kinematics, so theta_lab(theta_cm) and the recoil energy
   // both shift -- which is exactly why acceptance has to be measured per level, not assumed
   // from the elastic case.
   ExE.push_back(resEx);

   // ---- Recoil : p ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.0078250322); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg"
             << "   (recoil proton theta_lab ~ " << 0.5 * (180 - ThetaMaxCMS) << " - "
             << 0.5 * (180 - ThetaMinCMS) << " deg)" << std::endl;

   AtTPC2Body *TwoBody = new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE, ResEner,
                                        ThetaMinCMS, ThetaMaxCMS);
   primGen->AddGenerator(TwoBody);

   run->SetGenerator(primGen);
   run->SetStoreTraj(kTRUE);

   run->Init();

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);
   run->CreateGeometryFile("./data/geofile_C14_pp_full.root");

   timer.Stop();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << endl << endl;
}
