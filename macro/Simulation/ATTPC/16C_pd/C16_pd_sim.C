/// @file C16_pd_sim.C
/// @brief a1975 16C(p,d)15C simulation, matched to the experimental running conditions.
///
/// Purpose: the acceptance of the deuteron channel, measured per FINAL STATE. The (p,d)
/// analysis extracts four states of 15C -- the ground state, 0.740, 3.10 and 4.66 MeV -- and a
/// transfer reaction changes its two-body kinematics with every one of them, so theta_lab
/// (theta_cm) and the deuteron energy both move. An acceptance measured on one state does not
/// describe another; that is the whole reason for running this per level rather than once.
///
/// CONDITIONS (a1975, proton-target runs):
///     B = 2.85 T,  H2 at 300 torr,  drift length 1000 mm
///     beam 16C at 192 MeV = 12.0 MeV/u
///
/// !! THE GAS IS 300 torr, AND THE PRODUCTION RECONSTRUCTION DISAGREES !!
/// rho = 3.308e-5 g/cm3 (H2, 300 torr, 293 K), geometry ATTPC_H300torr_RT.root. But the a1975
/// data reconstruction uses ATTPC_H1bar_geomanager.root (8.27e-5) and fitUKF_a1975.C defaults to
/// gasDensity = 9.0e-5, i.e. 2.5 TIMES too much material. Energy loss scales with path length, so
/// that error is angle- and vertex-dependent rather than a constant offset. This is the same
/// class of mistake as the a1954 600-versus-300 torr error, which cost an entire production.
///
/// The 9.0e-5 figure is doubly wrong: even at 1 bar it is the 273 K density, where 293 K gives
/// 8.27e-5. Use the _RT media entries, which carry the room-temperature values.
///
/// BEAM ENERGY. 192 MeV is what the (p,d) analysis macros use. The (p,p') calibration of the
/// same data set prefers 195.5 MeV; the difference moves pz by 0.9 percent and is left as a
/// parameter rather than silently chosen here.
///
/// Beam momentum: the ion generator multiplies pz by A, so pz is GeV/c per nucleon.
///     m(16C) = 16.0147 u = 14917.60 MeV
///     KE = 192.0 MeV  ->  p = 2401.088 MeV/c  ->  pz = 2.401088/16 per nucleon
///
///   root -l 'C16_pd_sim.C(2000)'                      // ground state
///   root -l 'C16_pd_sim.C(8000,2.,178.,"TGeant4",-28.5,"./data/pd_ex1.root",0.740,1001)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. The discrepancy is
/// repo-wide and is kept so this sim stays comparable to its siblings -- but it means the beam
/// energy must be VERIFIED from MC truth, not assumed.
///
/// CM ANGULAR RANGE: the default is the full 2-178 deg, because an acceptance measurement has to
/// measure where the detector stops accepting. Truncating the generated range to the region that
/// is known to work builds the answer into the input.
void C16_pd_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 192.0)
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed: gRandom->GetSeed() returns TRandom3's internal state counter, not
   // the seed, so it cannot be used to verify that a parallel run was actually seeded.
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
   ATTPC->SetGeometryFileName("ATTPC_H300torr_RT.root"); // H2 at 300 torr, rho = 3.308e-5 g/cm3 (293 K)
   run->AddModule(ATTPC);

   // Field sign follows the a1954 lesson: generate in the DATA convention so that the sense of
   // rotation matches, otherwise anything inferring direction from curvature (AtSpyralPID) rejects
   // most of the simulated tracks and the acceptance measures the sign error instead of the
   // detector.
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
   Double_t Bmass = 16.0147;   // amu (repo convention -- see the note above)
   Double_t NomEnergy = Ebeam; // MeV

   // maxELoss CONTROLS THE REACTION RATE and is not cosmetic. AtTPCIonGenerator samples a
   // threshold uniformly in [0, maxELoss] and AtTpc::reactionOccursHere() fires only once the
   // accumulated loss exceeds it. Set it far above the real traversal loss and most events never
   // react; set it far below and every reaction happens near the entrance, so the vertex is no
   // longer uniform in z and the acceptance is measured on the wrong vertex distribution.
   //
   // MEASURED, not estimated: solving the elastic truth for the beam energy that gives Ex = 0
   // event by event showed 25.2 MeV of loss across the chamber when the gas was (wrongly) set to
   // 1 bar. Scaling by the density ratio 3.308/8.27 gives about 10 MeV at 300 torr. Re-verify the
   // same way after any change: the vertex must stay flat across the drift length, and the solved
   // beam energy should fall linearly from the entrance value.
   Double_t maxELoss = 10.0; // MeV

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 16C(p,d)15C   --------------------------------
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

   // ---- Target : p ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.0078250322); // uma
   ExE.push_back(0.0);

   // ---- Heavy residual : 15C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(15);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(15.0105993); // uma
   // The four states the (p,d) analysis extracts: 0 (g.s., 1/2+), 0.740 (5/2+), 3.10, 4.66.
   // Note that 3.10 and 4.66 lie above the neutron separation energy of 15C (1.218 MeV), so the
   // residual is unbound and will decay. That is not simulated here: the residual is transported
   // as a stable ion. It does not affect the DEUTERON kinematics, which is what the acceptance is
   // measured on, but it does mean this simulation must not be used to say anything about the
   // recoil or about coincidences with the neutron.
   ExE.push_back(resEx);

   // ---- Light ejectile : d, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(2);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.0141017781); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "16C(p,d)15C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "residual 15C excitation = " << resEx << " MeV" << std::endl;
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
   run->CreateGeometryFile("./data/geofile_C16_pd_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
