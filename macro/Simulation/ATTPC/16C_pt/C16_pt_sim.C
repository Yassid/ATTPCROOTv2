/// @file C16_pt_sim.C
/// @brief a1975 16C(p,t)14C simulation, matched to the proton-target running conditions.
///
/// Purpose: the acceptance of the TRITON channel, and -- more urgently -- ground truth for the
/// PID gate. 16C(p,t)14C has no measured signal in a1975: the excitation peak WALKS about
/// 11.5 MeV between theta_lab 8 and 40 deg, which is the signature of assigning triton mass to
/// particles that are not tritons. The gate was drawn by eye on a plane with no MC truth behind
/// it, and with an expected signal of only 1-10% of the (p,d) yield sitting under (p,p) and (p,d)
/// channels that are 10-100x more numerous, a gate that leaks 1% of the deuterons swamps it.
/// (d,t) settled the same question by drawing its gate on the SIMULATION's own plane and scoring
/// it against MC truth: 99.7% efficiency, 99.7% purity. That is what this exists to make possible.
///
/// CONDITIONS (a1975, proton-target runs 0106-0189 -- the same block as (p,d) and (p,p)):
///     B = 2.85 T,  H2 at 300 torr (rho 3.308e-5 g/cm3, 293 K),  beam 16C at 185 MeV = 11.6 MeV/u
///
/// BEAM ENERGY IS 185, NOT THE (p,d) SIM'S 192. working_point.sh exports A1975_PT_EBEAM = 185.0,
/// which is what the (p,t) analysis was done at. The (p,d) simulation macro defaults to 192 while
/// the (p,d) ANALYSIS uses 183.5-185; that inconsistency is documented in working_point.sh and is
/// not inherited here. Generate at the energy the analysis uses, or the acceptance describes a
/// different experiment.
///
/// THE REACTION IS TWO-NEUTRON PICKUP AND ITS Q VALUE IS POSITIVE:
///     Q = B(t) - S_2n(16C) = 8.482 - 5.469 = +3.013 MeV
/// against -2.03 for (p,d). 16C is neutron-rich, so taking two neutrons is energetically
/// downhill, and the tritons come out ENERGETIC: 4.7 to 29.5 MeV over theta_lab 3-37 deg for the
/// ground state. That matters for the fitter -- genfit applies NO stopping power below
/// beta*gamma = 0.05, which is 3.5 MeV for a triton, and unlike (d,t) (whose tritons reach down to
/// 0.9 MeV) this channel sits entirely above that threshold.
///
/// THE STATES. 14C has nothing between the 0+ ground state and 6.09 MeV, then a cluster --
/// 6.09 (1-), 6.59 (0+), 6.73 (3-), 7.01 (2+) -- that this resolution (~0.6 MeV FWHM) cannot
/// separate. All five are generated anyway: not to resolve them, but because Ex moves the
/// theta_lab(theta_cm) mapping, and generating the cluster lets the state dependence of the
/// acceptance be MEASURED rather than assumed. (d,t) found its acceptance state-independent to
/// 1%; whether that holds here is a result, not an input.
///
/// CM ANGULAR RANGE: the full 2-178 deg, and this is a real difference from (d,t). That macro
/// truncated at 70 because its high-KE branch has an ~840 mm cyclotron radius in a ~290 mm
/// chamber and is simply not measurable. Here the radius peaks at about 298 mm near theta_cm 80-90
/// and comes back down, so the backward branch IS potentially in acceptance -- and for a PURITY
/// study it is exactly the branch worth having, because those are the stiff tracks most easily
/// confused with something else. Truncating the generated range to the region known to work builds
/// the answer into the input.
///
/// AtTPC2Body ranges over the RESIDUAL's cm angle, and the analysis convention is
/// theta_cm = pi - theta3_cm ("the residual recoils opposite the ejectile"). theta_cm 2-70 maps to
/// triton theta_lab 3-36 deg; 70-178 folds BACK down to 0.8 deg on the high-KE branch, so a given
/// lab angle carries two solutions. That ambiguity is physics, not a bug, and it is one more thing
/// the PID gate has to survive.
///
///   root -l 'C16_pt_sim.C(2000)'                       // ground state
///   root -l 'C16_pt_sim.C(20000,2.,178.,"TGeant4",-28.5,"./data/pt_ex1.root",6.09,3001)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. The discrepancy is
/// repo-wide and is kept so this sim stays comparable to its siblings -- but it means the beam
/// energy must be VERIFIED from MC truth, not assumed.
void C16_pt_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 185.0)
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed: gRandom->GetSeed() returns TRandom3's internal state counter, not
   // the seed, so it cannot be used to verify that a parallel run was actually seeded.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar_pt.root";

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

   // -----   Two-body reaction 16C(p,t)14C   --------------------------------
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

   // ---- Heavy residual : 14C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(14);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(14.0032420); // uma
   // 14C states: the 0+ ground state, then 6.09 (1-), 6.59 (0+), 6.73 (3-), 7.01 (2+). Every one
   // of the excited states is ABOVE the alpha threshold of 14C, so the residual is unstable and
   // will decay. That is not simulated here -- it is transported as a stable ion. It does not
   // affect the TRITON kinematics, which is what the acceptance is measured on, but it means this
   // simulation must not be used for anything about the recoil or about coincidences with it.
   ExE.push_back(resEx);

   // ---- Light ejectile : t, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(3);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(3.0160492779); // uma -- ATOMIC 3H, matching the atomic masses used above
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "16C(p,t)14C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "residual 14C excitation = " << resEx << " MeV" << std::endl;
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
   run->CreateGeometryFile("./data/geofile_C16_pt_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
