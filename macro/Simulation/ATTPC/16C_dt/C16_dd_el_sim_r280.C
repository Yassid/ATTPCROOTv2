/// @file C16_dd_el_sim_r280.C
/// @brief r280 copy (2026-11-10): gas to r 280 mm (ATTPC_D300torr_r280.root, as C16_dp_sim_r280.C) for the 16C(d,d)
///        elastic ACCEPTANCE of the scale test (D2_UKF/trigger/dd). Otherwise identical to C16_dd_el_sim.C.
/// Original brief: a1975 d(16C,d) ELASTIC simulation, same deuterium conditions as C16_dt_sim.C (copied from it, 2026-11-01).
///
/// Purpose: a deuteron reference for the PID / Brho estimators (AtHelix3D circle, spiral fit, Spyral). The (d,t)
/// sims contain no deuterons, so the slow-deuteron behaviour of the estimators (the "bulge" the spiral fit puts in
/// the data deuteron band at Brho 0.22-0.30) could not be checked against truth. Generated UNIFORM in theta_cm (no
/// cross-section weighting): this is an estimator study, not an acceptance or a yield.
///
/// Kinematics: recoil deuteron KE = 72.7 MeV cos^2(phi_lab); AtTPC2Body ranges over the RESIDUAL (16C) cm angle, so
/// 2-70 deg gives phi_lab ~89-55 deg, Brho ~0.03-1.0 T m -- the whole data deuteron band.
/// All conditions, conventions and caveats (mass in amu, field sign, maxELoss) as in C16_dt_sim.C.
///
///   root -l -b -q 'C16_dd_el_sim.C(5000,2.,70.,"TGeant4",-28.5,"/mnt/f/a1975_C16_dd_sim/el_s5001_sim.root",0.,5001)'
void C16_dd_el_sim_r280(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim_dd_el.root",
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
   TString parFile = outFile; parFile.ReplaceAll("_sim.root", "_par.root"); // next to the output, not ./data

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
   ATTPC->SetGeometryFileName("ATTPC_D300torr_r280.root"); // r 280 mm, as C16_dp_sim_r280.C
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

   // ---- Heavy residual : 16C (elastic) ----
   Zp.push_back(6);
   Ap.push_back(16);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(16.0147); // uma
   // Every state above 1.218 MeV (Sn of 15C) is unbound and will decay. That is NOT simulated: the
   // residual is transported as a stable ion. It does not affect the TRITON kinematics, which is
   // what the acceptance is measured on, but this simulation must not be used for the recoil or
   // for any coincidence with the neutron.
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
   std::cout << "d(16C,d) elastic  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "residual 16C excitation = " << resEx << " MeV" << std::endl;
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
   TString geoOut = outFile; geoOut.ReplaceAll("_sim.root", "_geofile.root");
   run->CreateGeometryFile(geoOut);

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
