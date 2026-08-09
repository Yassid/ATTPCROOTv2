/// @file C16_pp_a1975_sim.C
/// @brief a1975 16C(p,p')16C simulation, matched to the experimental running conditions.
///
/// Purpose: the ELASTIC channel, which is the reference the other channels are calibrated
/// against. Two questions it answers that the data alone cannot:
///
///   1. How much of the residual +0.20 MeV zero-point offset in the (p,p) spectrum is simply the
///      beam losing energy before it reacts? pp/ebeam_pp.C found that no beam energy between 186
///      and 202 MeV puts the ground state at Ex = 0 -- it runs +0.16 to +0.23 across that whole
///      range -- while the 1.766 MeV level spacing is matched at 195.3. One parameter cannot
///      satisfy both conditions, and beam energy loss is the obvious candidate for the residual,
///      because it shifts the reconstructed Ex without being a property of the nominal energy.
///   2. What the excitation-energy resolution is, and how it varies with angle.
///
/// CONDITIONS (a1975, proton-target runs):
///     B = 2.85 T,  H2 at 1 bar,  drift length 1000 mm
///     beam 16C at 192 MeV = 12.0 MeV/u
///
/// THE GAS. The workspace documents H2 at 1 bar and the fitters are run with
/// gasDensity = 9.0e-5 g/cm3. That number is the 273 K density; at the 293 K the gas was
/// actually at, 1 bar of H2 is 8.27e-5 g/cm3, which is what the geometry ATTPC_H1bar.root and
/// the medium H_1bar carry and what is used here. The 9.0e-5 in the reconstruction is therefore
/// about 9 percent too dense. This is the same trap that was found and fixed in the a1954
/// analysis, where media.geo carried the 273 K value for 300 torr; it is recorded here because
/// energy loss scales with path length, so an error in it is angle- and vertex-dependent rather
/// than a constant offset.
///
/// Separately, the a1975 reconstruction loads parameters/ATTPC.a1954.par, which declares
/// GasPressure 600 torr and Density 0.197 mg/cm3 -- neither of which is H2 at 1 bar. Those
/// entries drive the digitisation, not the transport, but the inconsistency is real and worth
/// resolving before these acceptances are used for anything absolute.
///
/// BEAM ENERGY. 192 MeV is what the (p,d) analysis macros use. The (p,p') calibration of the
/// same data set prefers 195.5 MeV; the difference moves pz by 0.9 percent and is left as a
/// parameter rather than silently chosen here.
///
/// Beam momentum: the ion generator multiplies pz by A, so pz is GeV/c per nucleon.
///     m(16C) = 16.0147 u = 14917.60 MeV
///     KE = 192.0 MeV  ->  p = 2401.088 MeV/c  ->  pz = 2.401088/16 per nucleon
///
///   root -l 'C16_pp_a1975_sim.C(2000)'                      // elastic
///   root -l 'C16_pp_a1975_sim.C(8000,2.,178.,"TGeant4",-28.5,"./data/pp_2plus.root",1.766,1001)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. The discrepancy is
/// repo-wide and is kept so this sim stays comparable to its siblings -- but it means the beam
/// energy must be VERIFIED from MC truth, not assumed.
///
/// CM ANGULAR RANGE: the default is the full 2-178 deg, because an acceptance measurement has to
/// measure where the detector stops accepting. Truncating the generated range to the region that
/// is known to work builds the answer into the input.
void C16_pp_a1975_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
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
   ATTPC->SetGeometryFileName("ATTPC_H1bar.root"); // H2 at 1 bar, rho = 8.27e-5 g/cm3 (293 K)
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
   // 16C at 12 MeV/u through 1 m of H2 at 1 bar loses roughly 28 MeV: the a1954 14C beam lost
   // ~12 MeV through 1 m at 300 torr, and this gas is 2.3x denser at a similar velocity.
   // VERIFY THIS on the first run by checking that the truth vertex z is flat across the drift
   // length -- a slope means the estimate is wrong.
   Double_t maxELoss = 28.0; // MeV

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

   // ---- Heavy residual : 16C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(16);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(16.0147); // uma
   // 0 = elastic; 1.766 = the first 2+ of 16C, which is the level the beam-energy calibration
   // uses as its ruler.
   ExE.push_back(resEx);

   // ---- Light ejectile : p, the particle that is detected ----
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
   std::cout << "16C(p,p')16C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
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
   run->CreateGeometryFile("./data/geofile_C16_pp_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
