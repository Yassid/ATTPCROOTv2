/// @file C15_dd_sim.C
/// @brief a2091 15C(d,d') simulation for the ACCEPTANCE and the RESPONSE FUNCTION of the
///        levels_15C analysis. Adapted from 16C_dt/C16_dd_el_sim.C (a1975), same gas and field.
///
/// WHY THIS EXISTS.  Every open item in macro/Unpack_HDF5/C15p/levels_15C/LEVELS_15C.md needs the
/// same missing ingredient:
///   - the continuum cannot be MODELLED, only absorbed into an empirical shape, because the
///     efficiency as a function of Ex is unknown (the kinematic acceptance is computed exactly in
///     binacc.h and is ~1 inside a theta_CM bin, so what is missing is the RECONSTRUCTION
///     efficiency, not the geometry);
///   - the response function is a bifurcated Gaussian with exponential tails fitted to the data
///     itself, so the quoted FWHMs are observed widths and the 5/2+ yield is lineshape-limited;
///   - theta_CM bin 0 is degenerate between the background and BW1 and cannot be used for yields.
///
/// The residual 15C is transported as a STABLE ion even above Sn = 1.218 MeV, exactly as in the
/// a1975 macro. That is correct for the deuteron kinematics, which is all the acceptance is
/// measured on, but this simulation must not be used for the neutron or any coincidence.
///
/// Ebeam defaults to 191.25 MeV, the value the levels_15C analysis uses for ATTPCROOT, so the
/// simulated and measured spectra share one energy scale. Pass 195.0 for the Spyral comparison.
///
///   root -l -b -q 'C15_dd_sim.C(5000, 2., 70., "TGeant4", -28.5, "out_sim.root", 0.0, 5001, 191.25)'
///                                                                               ^Ex of the residual
void C15_dd_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim_C15dd.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 191.25)
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed: gRandom->GetSeed() returns TRandom3's internal state counter, not
   // the seed, so it cannot be used to verify that a parallel run was actually seeded.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   // Derive the par and geometry files from outFile so a parallel job writes them NEXT TO its
   // own output, not into a ./data/ that may not exist (FairParRootFileIo segfaults rather than
   // reporting a missing directory).
   TString parFile = outFile; parFile.ReplaceAll("_sim.root", "_par.root");
   TString geoFile = outFile; geoFile.ReplaceAll("_sim.root", "_geo.root");

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
   Int_t a = 15; // mass number
   Int_t q = 0;  // charge state
   Int_t m = 1;  // multiplicity
   const Double_t u = 931.49401, mBeamMeV = 15.0105993 * u;
   Double_t pTot = std::sqrt(Ebeam * Ebeam + 2 * Ebeam * mBeamMeV) / 1000.0; // GeV/c
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = pTot / a; // GeV/c per nucleon
   Double_t BExcEner = 0.0;
   Double_t Bmass = 15.0105993; // amu (repo convention -- see the note above)
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

   // -----   Two-body reaction 15C(d,d')15C*   --------------------------------
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
   Mass.push_back(15.0105993); // uma
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

   // ---- Heavy residual : 15C, excited to resEx ----
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
   std::cout << "d(15C,d')  Ebeam = " << Ebeam << " MeV (" << Ebeam / 15 << " MeV/u),  Ex(res) = " << resEx << " MeV,  p = " << pTot * 1000
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
   run->CreateGeometryFile(geoFile.Data());

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
