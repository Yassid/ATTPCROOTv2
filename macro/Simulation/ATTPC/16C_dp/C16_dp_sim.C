/// @file C16_dp_sim.C
/// @brief a1975 16C(d,p)17C simulation, matched to the deuterium running conditions.
///
/// Purpose: the acceptance AND the resolution floor of the proton channel, per final state.
/// Written 2026-09-23 to answer a specific question the data could not settle on its own --
/// the BOUND region of the (d,p) Ex spectrum shows, under the adopted selection:
///     * a low-Ex peak that WALKS +0.30 MeV from one end of the chamber to the other, which is
///       the size of the entire 17C bound level scheme (g.s. / 0.217 / 0.331), and
///     * a resolution collapse at v_z 650-750 mm (sigma 0.25 -> 0.93 MeV) that empties the
///       bound window there while the total yield stays flat.
/// Neither can be called geometry or reconstruction without a simulation carrying truth through
/// the SAME chain. That is what this is for; see D2_UKF/exsel_dp.C for the selection it must be
/// pushed through.
///
/// THIS IS A DIRECT COPY OF C16_dt_sim.C with the reaction products swapped. Everything about the
/// running conditions is deliberately identical, so the two channels' acceptances are comparable:
///     B = 2.85 T,  D2 at 300 torr,  drift velocity 1.10424 cm/us, drift length 971.7312 mm
/// Read that file's header for why each choice is what it is -- it is not repeated here.
///
/// BEAM ENERGY 184.25 MeV, chosen by Yassid 2026-09-23. The (d,p) channel has carried THREE
/// values on the same 47 runs -- 184.25 (working_point.sh, measured on these runs and adopted for
/// (d,t)), 188.33 (the old Spyral-era C17_dp_fits analysis, 11.77 MeV/u) and 192.0 (an old macro
/// default). 184.25 is the one the data selection under study uses, so it is the only one that
/// makes sim and data comparable. It is NOT independently settled -- see project_a1975_dp_analysis.
///
/// MASS CONVENTION: ATOMIC masses throughout, which is what the sibling sim macros pass and what
/// AtTPCIonGenerator/AtTPC2Body expect. This is NOT the nuclear-light/atomic-heavy convention the
/// ANALYSIS uses when it inverts kinematics (see reference_mass_convention_electrons) -- do not
/// copy numbers between the two without checking which convention they are in.
///
/// SIM z HANDEDNESS: the simulation reverses drift-z against the experiment, so fitting this
/// output needs bFieldSign = +1 where the data uses -1, and the reconstructed polar angle comes
/// back as 180 - true. That is a property of digitisation; do not "fix" it here.
///
/// CM ANGULAR RANGE: default is the full 2-178 deg. An acceptance has to measure where the
/// detector stops accepting; generating only where it is known to work builds the answer in.
/// NOTE the convention -- AtTPC2Body ranges over the RESIDUAL's cm angle, so the ejectile angles
/// that come out are pi - theta_cm. Verify against truth before trusting any narrowed range.
///
///   root -l 'C16_dp_sim.C(20000)'                                   // ground state, full range
///   root -l 'C16_dp_sim.C(20000,2.,178.,"TGeant4",-28.5,"./data/dp_gs.root",0.0,4001,184.25)'
void C16_dp_sim(Int_t nEvents = 20000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim_dp.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 184.25)
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar_dp.root";

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
   ATTPC->SetGeometryFileName("ATTPC_D300torr_v2.root");
   run->AddModule(ATTPC);

   // Field sign follows the data convention so the sense of rotation matches; otherwise anything
   // inferring direction from curvature rejects most simulated tracks and the acceptance measures
   // the sign error instead of the detector.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 16C   ---------------------------------------------------
   Int_t z = 6;  // atomic number
   Int_t a = 16; // mass number
   Int_t q = 0;  // charge state
   Int_t m = 1;  // multiplicity
   const Double_t u = 931.49401, mBeamMeV = 16.0147 * u;
   Double_t pTot = std::sqrt(Ebeam * Ebeam + 2 * Ebeam * mBeamMeV) / 1000.0; // GeV/c
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = pTot / a; // GeV/c per nucleon
   Double_t BExcEner = 0.0;
   Double_t Bmass = 16.0147;   // amu (repo convention -- the mass argument is documented as GeV,
                               // the repo passes amu; verify the beam energy from MC truth)
   Double_t NomEnergy = Ebeam; // MeV

   // maxELoss CONTROLS THE REACTION RATE and is not cosmetic -- it sets where along z the reaction
   // fires, hence the VERTEX DISTRIBUTION, which is exactly what this study is about. 10 MeV is the
   // (p,d)/(d,t) value and carries over (same number density, same stopping).
   // *** VERIFY on the output that the vertex stays flat along the drift length before reading any
   // v_z-dependent result out of this simulation. A sloped generated vertex distribution would
   // imitate the very effect being hunted. ***
   Double_t maxELoss = 10.0; // MeV

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 16C(d,p)17C   --------------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, heavy residual, light ejectile -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 16C ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(16.0147); // uma
   ExE.push_back(BExcEner);

   // ---- Target : d ----
   Zp.push_back(1);
   Ap.push_back(2);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.0141017781); // uma
   ExE.push_back(0.0);

   // ---- Heavy residual : 17C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(17);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(17.0225864); // uma
   // Sn(17C) = 0.729 MeV, so every state above that is unbound and will decay. That is NOT
   // simulated: the residual is transported as a stable ion. It does not affect the PROTON
   // kinematics, which is what the acceptance is measured on, but this simulation must not be used
   // for the recoil or for any coincidence with the neutron.
   ExE.push_back(resEx);

   // ---- Light ejectile : p, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.00782503207); // uma (atomic H -- see the mass-convention note above)
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "16C(d,p)17C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "residual 17C excitation = " << resEx << " MeV" << std::endl;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg" << std::endl;

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
   run->CreateGeometryFile("./data/geofile_C16_dp_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
