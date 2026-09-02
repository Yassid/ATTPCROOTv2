/// @file C17_dp_sim.C
/// @brief 17C(d,p)18C simulation for the SOLARIS + AT-TPC 17C proposal (C17p_FRIB_Proposal).
///
/// WHY THIS EXISTS. The proposal's primary goal is M_n/M_p for the 1/2+ (217 keV) and 5/2+
/// (332 keV) states of 17C, from 17C(p,p') and 17C(d,d') on hydrogen and deuterium gas. It also
/// notes that the deuterium run gives (d,p), (p,d), (d,3He) etc. "for free", as additional
/// structure information on the neutron-rich carbons. This macro is the (d,p) arm of that: what
/// 17C(d,p)18C looks like in the AT-TPC during the deuterium day.
///
/// It is the 14C(d,p)15C simulation with the beam, the target-residual pair and the level scheme
/// changed, and NOTHING else -- same gas, same field, same drift, same beam hole, same PSA and
/// HDBSCAN settings, same genfit configuration. That is deliberate: 14C(d,p) is a measured,
/// debugged reference (see macro/Simulation/ATTPC/14C_dp/RESULTS.md), so every difference between
/// the two is the reaction and not a setting.
///
/// CONDITIONS (from the proposal, ReA6):
///     beam    17C at 8.37 MeV/u = 142.29 MeV, 940 pps
///     field   B = 2.85 T, the nominal SOLARIS + AT-TPC field and the one the 14C(p,p')
///             measurement of Ref. [24] (Ayyad et al., Front. Phys. 13, 1539148 (2025)) ran at
///     gas     D2 at 300 torr, 293 K -- the deuterium counterpart of the 300 torr H2 that the
///             same 14C measurement used
///     drift   1000 mm
///
/// ANGULAR DISTRIBUTION IS FLAT. AtTPC2Body samples theta_cm flat in COS(theta_cm), i.e.
/// isotropic in the centre of mass (AtTPC2Body.cxx:169). There is no DWBA angular distribution
/// here and none is wanted yet: with a flat generator the per-bin acceptance is a clean ratio and
/// the Ex resolution is measured per angular slice, so both are independent of the true shape.
/// When the DWBA calculations exist, they weight these results; they do not require regenerating.
///
/// LEVELS. 18C has S_n = 4184 keV, and the ENSDF adopted level scheme has exactly FOUR bound
/// states below it:
///     0      0+          ground state
///     1588   2+          T_1/2 = 15.5 ps, B(E2) = 0.000364 e2b2 -- the state the proposal's
///                        B(E2) discussion is about
///     2515   (2+)        T_1/2 < 3.2 ps
///     3972   (2,3)+
/// All four are simulated. The next thing up is the neutron continuum, so this is the whole bound
/// spectrum -- unlike the 14C(d,p) study, which had to simulate unbound 15C levels to get enough
/// lever arm on the resolution. 927 keV between the 1588 and the 2515 is the separation that
/// matters, and the proposal quotes ~300 keV as the achievable AT-TPC resolution.
///
/// Q VALUE. Q(17C(d,p)18C) = +1.959 MeV, from the AME2020 mass excesses
///     17C 21031.880   d 13135.723   18C 24919.266   p 7288.971  (keV)
/// This is POSITIVE, where 14C(d,p)15C is Q = -1.007 MeV. So at the same theta_cm the 17C protons
/// are ~3 MeV faster than their 14C counterparts, on a beam that is 3.1 MeV/u slower. The two
/// effects push the proton energy in opposite directions and the net is not obvious a priori --
/// which is one of the things this simulation is for.
///
/// Beam momentum: the ion generator multiplies the pz argument by A, so pz is GeV/c per nucleon
/// and the TOTAL momentum is what sets the beam energy.
///     m(17C) = 17.02257865 u = 15.8564298 GeV
///     KE = 142.29 MeV  ->  p = 2.1290066 GeV/c  ->  pz = 2.1290066/17 per nucleon
/// NomEnergy does NOT set the beam energy -- only pz does; AtVertexPropagator::GetBeamNomE() has
/// no consumer anywhere in the tree. But it is NOT harmless either: the same constructor argument
/// doubles as the default for maxELoss, which controls how often a reaction is generated at all
/// (see the maxELoss comment below). Pass maxELoss explicitly.
///
///   root -l 'C17_dp_sim.C(2000)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. That discrepancy is
/// repo-wide and is kept here so this sim stays comparable to its siblings -- but it means the
/// beam energy MUST be verified from MC truth rather than assumed. Run
/// 14C_pp/check_beam_C14.C on the output; it is beam-agnostic.

/// CM ANGULAR RANGE: 2-178 by default, i.e. everything. As in 14C(d,p), the (d,p) proton does NOT
/// follow theta_lab = (180 - theta_cm)/2: the transfer peak sits at small theta_cm, where the
/// proton comes out BACKWARD in the lab. Do not preselect -- the acceptance is measured per bin.
void C17_dp_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                Double_t resEx = 0.0, UInt_t seed = 0,
                // TRANSPORT GAS. D2 at 300 torr, 293 K, rho = 6.5643e-5 g/cm3. This is matched
                // EXACTLY to the digitisation par (ATTPC.C17dp_D300torr_b285.par, Density
                // 0.065643 mg/cm3) so that transport and digitisation model the same gas. Getting
                // that wrong is never inert here: the transport geometry sets the energy loss that
                // goes into the TRUTH.
                TString geoFile = "ATTPC_D300torr_v2.root")
{
   // RNG seed. seed = 0 keeps ROOT's default (time-based) behaviour; pass a distinct value per
   // parallel job, or byte-identical events get counted as added statistics.
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
   ATTPC->SetGeometryFileName(geoFile); // default ATTPC_D300torr_v2.root, rho = 6.5643e-5 g/cm3
   std::cout << "\033[1;33m[C17_dp_sim] transport gas geometry: " << geoFile << "\033[0m" << std::endl;
   run->AddModule(ATTPC);

   // -----   Magnetic field : 2.85 T   --------------------------------------
   // SIGN: generate in the DATA convention, which is NEGATIVE. AtSpyralPID::Direction() infers
   // forward/backward purely from the sense of rotation (it never looks at z) while `polar` comes
   // from d(rho)/dz, so generating at +28.5 kG makes its consistency check reject ~82 % of sim
   // tracks. This was measured on a1954 and it is not specific to that beam.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 17C at 142.29 MeV (8.37 MeV/u)   ------------------------
   Int_t z = 6;  // Atomic number
   Int_t a = 17; // Mass number
   Int_t q = 0;  // Charge state
   Int_t m = 1;  // Multiplicity
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = 2.1290066 / a; // GeV/c per nucleon -> 142.29 MeV total KE
   Double_t BExcEner = 0.0;
   Double_t Bmass = 17.02257865; // amu, AME2020 (repo convention -- see the note above)
   Double_t NomEnergy = 142.29;  // MeV

   // maxELoss CONTROLS THE REACTION RATE -- it is not cosmetic.
   // AtTPCIonGenerator does  fMaxEnLoss = (eLoss < 0 ? ener : eLoss)  and then samples
   //     Er = gRandom->Uniform(0, fMaxEnLoss);  SetRndELoss(Er)
   // while AtTpc::reactionOccursHere() fires only once  fELossAcc > GetRndELoss().
   // Leaving eLoss at its -1 default makes the threshold uniform in [0, 142.29], which exceeds the
   // available energy loss ~90 % of the time and generates no reaction.
   // MEASURED with CATIMA (AtELossCATIMA, D2 at 6.5643e-5 g/cm3): a 17C beam at 142.29 MeV loses
   // 14.80 MeV across the full 1000 mm (142.29 -> 127.49 MeV, i.e. 8.37 -> 7.50 MeV/u), at
   // 14.18 MeV/m entering. Its full range is 5.85 m, so the beam crosses the chamber comfortably.
   // The same calculation gives 11.21 MeV for the 14C(d,p) case, against the 12.0 that macro uses,
   // so the two are on the same footing.
   Double_t maxELoss = 15.0; // MeV, ~ the 17C energy loss across the full drift length

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 17C(d,p)18C   --------------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, scattered, recoil -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 17C ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(17.02257865); // uma
   ExE.push_back(BExcEner);

   // ---- Target : d ----
   Zp.push_back(1);
   Ap.push_back(2);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.0141017778); // uma
   ExE.push_back(0.0);

   // ---- Scattered : 18C (the residual, and the one that carries the excitation) ----
   Zp.push_back(6);
   Ap.push_back(18);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(18.02675193); // uma, AME2020
   // Excitation of the RESIDUAL 18C. 0 = the 0+ ground state; 1.588 = the 2+; 2.515 = the (2+);
   // 3.972 = the (2,3)+. All four are BOUND (S_n = 4.184 MeV).
   // The excitation changes the two-body kinematics, so theta_lab(theta_cm) and the proton energy
   // both shift -- which is exactly why acceptance has to be measured per level, not assumed.
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
             << "   (flat in cos(theta_cm); the (d,p) proton does NOT map as (180-theta_cm)/2)" << std::endl;

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
   run->CreateGeometryFile("./data/geofile_C17_dp_full.root");

   timer.Stop();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << endl << endl;
}
