/// @file C16_dt_sim.C
/// @brief a1975 16C(d,t)15C simulation, matched to the deuterium running conditions.
///
/// Purpose: the acceptance of the triton channel, per FINAL STATE. The (d,t) analysis extracts
/// five states of 15C -- g.s., 0.740, 3.103, 4.657 and a feature near 6.6 -- and a transfer
/// reaction changes its two-body kinematics with every one of them, so theta_lab(theta_cm) and
/// the triton energy both move. An acceptance measured on one state does not describe another.
///
/// CONDITIONS (a1975, deuterium-target runs). These are NOT the (p,d) conditions:
///     B = 2.85 T,  D2 at 300 torr (rho 6.61e-5 g/cm3),  beam 16C at 184.17 MeV = 11.5 MeV/u
///     drift velocity 1.10424 cm/us, drift length 971.7312 mm, TB entrance 560
/// The gas is twice as dense as the (p,d) hydrogen target by MASS, but the stopping power is
/// nearly the same: D2 and H2 at the same pressure and temperature have the same molecular number
/// density, hence the same electron density, and dE/dx follows the electrons. So maxELoss carries
/// over from the (p,d) macro unchanged; it is the MASS density that doubles, and that only enters
/// through the geometry medium, which is already TargetD2_300.
///
/// WHY THIS CHANNEL IS HARD, and what the acceptance is expected to show: a 1 MeV triton has a
/// 144 mm range and an 81 mm cyclotron radius at 2.85 T, so it stops before completing an arc.
/// The (d,t) locus is DOUBLE-VALUED below the ~56 deg turnover and the data populates the LOW-KE
/// branch (Brho 0.22-0.58); the high branch would need a ~840 mm radius in a ~290 mm chamber and
/// is not measurable. A correct acceptance must therefore fall off hard at forward theta_cm,
/// where the triton is slowest -- that is the shape to check this simulation against.
///
/// Beam momentum: the ion generator multiplies pz by A, so pz is GeV/c per nucleon.
///     m(16C) = 16.0147 u = 14917.60 MeV,  KE = 184.17 MeV  ->  p = 2350.9 MeV/c
///
///   root -l 'C16_dt_sim.C(2000)'                        // ground state
///   root -l 'C16_dt_sim.C(8000,2.,178.,"TGeant4",-28.5,"./data/dt_ex1.root",0.740,1001)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. The discrepancy is
/// repo-wide and is kept so this sim stays comparable to its siblings -- but it means the beam
/// energy must be VERIFIED from MC truth, not assumed. dt_truth_check.C does exactly that.
///
/// CM ANGULAR RANGE: the default is the full 2-178 deg. An acceptance has to measure where the
/// detector stops accepting; generating only where it is known to work builds the answer in.
void C16_dt_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim_dt.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 184.17)
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed: gRandom->GetSeed() returns TRandom3's internal state counter, not
   // the seed, so it cannot be used to verify that a parallel run was actually seeded.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar_dt.root";

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
   // TargetD2_300 in media.geo: 6.5643e-5 g/cm3, i.e. D2 at 300 torr and about 295 K. The data
   // parameter file declares 6.61e-5 (293 K). The 0.7 % difference is left alone deliberately --
   // changing media.geo is dangerous (appending past Ar90CF4_250mbar hangs FairGeoMedia) and a
   // 0.7 % density error is far below the systematics this acceptance carries.
   ATTPC->SetGeometryFileName("ATTPC_D300torr_v2.root");
   run->AddModule(ATTPC);

   // Field sign follows the a1954 lesson: generate in the DATA convention so the sense of rotation
   // matches, otherwise anything inferring direction from curvature rejects most simulated tracks
   // and the acceptance measures the sign error instead of the detector.
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

   // maxELoss CONTROLS THE REACTION RATE and is not cosmetic. AtTPCIonGenerator samples a threshold
   // uniformly in [0, maxELoss] and AtTpc::reactionOccursHere() fires only once the accumulated loss
   // exceeds it. Too high and most events never react; too low and every reaction happens near the
   // entrance, so the vertex stops being uniform in z and the acceptance is measured on the wrong
   // vertex distribution. 10 MeV is the (p,d) value and carries over: same number density, same
   // stopping. VERIFY on the output that the vertex stays flat along the drift length.
   Double_t maxELoss = 10.0; // MeV

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 16C(d,t)15C   --------------------------------
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

   // ---- Heavy residual : 15C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(15);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(15.0105993); // uma
   // Every state above 1.218 MeV (Sn of 15C) is unbound and will decay. That is NOT simulated: the
   // residual is transported as a stable ion. It does not affect the TRITON kinematics, which is
   // what the acceptance is measured on, but this simulation must not be used for the recoil or
   // for any coincidence with the neutron.
   ExE.push_back(resEx);

   // ---- Light ejectile : t, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(3);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(3.0160492779); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "16C(d,t)15C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
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
   run->CreateGeometryFile("./data/geofile_C16_dt_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
