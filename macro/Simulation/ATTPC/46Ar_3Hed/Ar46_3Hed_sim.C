/// @file Ar46_3Hed_sim.C
/// @brief 46Ar(3He,d)47K at 13 MeV/u in the AT-TPC, active-target mode.
///
/// Proposal (ar46hd_attpc_sub.pdf) asks for the single-proton strength of 47K: the 1/2+ ground
/// state (1s1/2), the 0.36 MeV 3/2+ (0d3/2) and the 2.02 MeV 7/2- (0f7/2), over 15 < theta_cm <
/// 80 deg, aiming at ~350 keV FWHM in Q value.
///
/// GAS-CELL MODE WAS DROPPED. The proposal puts 300 torr of 3He in a 1-inch PPTA cell on the beam
/// axis and fills the tracking volume with isobutane. Here the whole drift volume IS the target:
/// 3He + 5% CO2 (molar) at 300 torr, 293.15 K, rho = 8.3128e-5 g/cm3, geometry
/// ATTPC_He3CO2_300torr.root. Over the 100 cm drift that is 4.70 mg/cm2 of 3He = 9.4e20 at/cm2,
/// i.e. the target thickness the proposal quotes for the cell, without the cell wall or the
/// isobutane -- so the luminosity estimate carries over and the vertex is measured, not assumed.
///
/// MEDIA ENTRY. He3CO2_300torr sits in media.geo ABOVE Ar90CF4_250mbar, deliberately. The last
/// entry in that file is malformed (it is below the "additional parameters" caution line but
/// carries no extra parameter line), which is harmless only while it is last: anything appended
/// after it -- media entries or even bare // comment lines -- desyncs FairGeoMedia, which then
/// hangs at "Read media" and allocates until the machine dies. Do not move this block down, and
/// do not append to media.geo.
///
/// CM ANGLE CONVENTION. AtTPC2Body samples fThetaCms for the HEAVY product (index 2) and gives
/// the light ejectile pi - theta (AtTPC2Body.cxx:209). That makes fThetaCms identical to the
/// proposal's DWBA theta_cm, so the arguments below are the proposal's numbers unmodified:
///
///     theta_cm    theta_lab(d)    T_d       helix diameter
///        15 deg     131.4 deg     5.1 MeV      24 cm
///        30 deg     104.3 deg    11.5 MeV      47 cm
///        50 deg      82.7 deg    25.6 MeV      72 cm
///        80 deg      59.8 deg    55.3 MeV      93 cm
///
/// which reproduces the proposal's own "130 (15) to 60 (80)" line. Q(gs) = +7.776 MeV. The
/// drift volume is only 50 cm across, so the forward-CM deuterons curl inside it while the
/// backward ones leave radially -- the acceptance turn-over is physics of the setup, so the
/// default generated range is the full 15-80 and the cut is made afterwards, not here.
///
/// Beam momentum: the ion generator multiplies pz by A, so pz is GeV/c per nucleon.
///     m(46Ar) = 45.9680367 u,  KE = 598 MeV (13 MeV/u)  ->  p = 7180.4 MeV/c
///
///   root -l 'Ar46_3Hed_sim.C(2000)'                                            // ground state
///   root -l 'Ar46_3Hed_sim.C(8000,15.,80.,"TGeant4",-28.5,"./data/d_360.root",0.360,1001)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. The discrepancy is
/// repo-wide and is kept so this sim stays comparable to its siblings -- but it means the beam
/// energy must be VERIFIED from MC truth, not assumed. Atomic masses are used throughout: 18+2
/// electrons on the entrance channel and 19+1 on the exit channel cancel exactly, so Q is right.
void Ar46_3Hed_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 15.0, Double_t thetaMaxCM = 80.0,
                   TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5,
                   TString outFile = "./data/attpcsim.root", Double_t resEx = 0.0, UInt_t seed = 0,
                   Double_t Ebeam = 598.0,
                   // FORWARD TELESCOPE, OFF BY DEFAULT. Adding a module changes the geometry, so
                   // switching it on makes new sims differ from every sim already on disk. It is
                   // opt-in for exactly that reason: an unqualified rerun must reproduce what is
                   // there. Pass kTRUE to place the two DSSDs behind the cathode -- only useful
                   // for the REVERSED detector, where the beam leaves through the cathode.
                   Bool_t withTelescope = kFALSE,
                   TString telescopeGeo = "Ar46_telescope_v1.0.root")
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
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
   ATTPC->SetGeometryFileName("ATTPC_He3CO2_300torr.root"); // 3He + 5% CO2, 300 torr, 8.3128e-5 g/cm3
   run->AddModule(ATTPC);

   // The dE-E telescope is a SEPARATE active module producing AtSiPoint, the same arrangement the
   // HELIOS macros use for AtSiArray. Its geometry (geometry/Ar46_telescope.C) places the two
   // DSSDs at z = 105 and 106 cm, i.e. past the 100 cm drift volume, so it cannot overlap the TPC.
   if (withTelescope) {
      FairDetector *tel = new AtSiArray("Ar46Telescope", kTRUE);
      tel->SetGeometryFileName(telescopeGeo);
      run->AddModule(tel);
      std::cout << "  FORWARD TELESCOPE ON: " << telescopeGeo << " (AtSiPoint branch will be written)\n";
   }

   // Field sign follows the a1954/a1975 lesson: generate in the DATA convention so the sense of
   // rotation matches, otherwise anything inferring direction from curvature rejects most of the
   // simulated tracks and the acceptance measures the sign error instead of the detector.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 46Ar   --------------------------------------------------
   Int_t z = 18; // atomic number
   Int_t a = 46; // mass number
   Int_t q = 0;  // charge state
   Int_t m = 1;  // multiplicity
   const Double_t u = 931.49401, mBeamMeV = 45.9680367 * u;
   Double_t pTot = std::sqrt(Ebeam * Ebeam + 2 * Ebeam * mBeamMeV) / 1000.0; // GeV/c
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = pTot / a; // GeV/c per nucleon
   Double_t BExcEner = 0.0;
   Double_t Bmass = 45.9680367; // amu (repo convention -- see the note above)
   Double_t NomEnergy = Ebeam;  // MeV

   // maxELoss CONTROLS THE REACTION RATE and is not cosmetic. AtTPCIonGenerator samples a
   // threshold uniformly in [0, maxELoss] and AtTpc::reactionOccursHere() fires only once the
   // accumulated loss exceeds it. Set it far above the real traversal loss and most events never
   // react; set it far below and every reaction happens near the entrance, so the vertex is no
   // longer uniform in z.
   //
   // MEASURED, not estimated. Running this macro with maxELoss = 1e5 (no reaction ever fires) and
   // summing fELoss of track 0 over the full drift length gives, per beam traversal:
   //
   //     deposit  fELossAcc          95.7 MeV   (0.957 MeV/cm)   <-- what the threshold sees
   //     true beam kinetic energy    98.6 MeV   (13.01 -> 10.86 MeV/u, from the point momenta)
   //     Bethe estimate                122 MeV  (no density/shell corrections, so it runs high)
   //
   // WATCH THE DENOMINATOR when you re-measure: AtTPCIonGenerator alternates beam events and
   // reaction events (AtTPCIonGenerator.cxx:166 adds the beam track only on beam events), so half
   // the tree entries carry no beam track at all. Averaging the deposit over ALL entries reads
   // 47.8 MeV instead of 95.7 and would send you looking for a factor of two that is not there.
   //
   // 96 MeV therefore puts the threshold at the far end of the chamber: essentially every
   // reaction event fires, with the vertex flat across the drift length.
   //
   // The kinematics are NOT affected by this knob: startReactionEvent() hands the generator
   // gMC->Etot() (AtTpc.cxx:257), the true transported energy, so the residual beam energy at the
   // vertex is right whatever maxELoss is. maxELoss only places the vertex.
   Double_t maxELoss = 96.0; // MeV

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 46Ar(3He,d)47K   -----------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, heavy residual, light ejectile -- must be 4
   Double_t ResEner = 0.0; // unused in active-target mode (taken from the vertex propagator)

   // ---- Beam : 46Ar ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(45.9680367); // uma
   ExE.push_back(BExcEner);

   // ---- Target : 3He ----
   Zp.push_back(2);
   Ap.push_back(3);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(3.01602932); // uma
   ExE.push_back(0.0);

   // ---- Heavy residual : 47K, left in the state under study ----
   Zp.push_back(19);
   Ap.push_back(47);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(46.9616160); // uma
   // 0 = 1/2+ g.s. (1s1/2); 0.360 = 3/2+ (0d3/2); 2.020 = 7/2- (0f7/2)
   ExE.push_back(resEx);

   // ---- Light ejectile : d, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(2);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.01410178); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "46Ar(3He,d)47K  Ebeam = " << Ebeam << " MeV (" << Ebeam / 46 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "47K excitation = " << resEx << " MeV" << std::endl;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg" << std::endl;

   AtTPC2Body *TwoBody =
      new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE, ResEner, ThetaMinCMS, ThetaMaxCMS);
   primGen->AddGenerator(TwoBody);

   run->SetGenerator(primGen);

   // Trajectory storage makes the output very large; for display only, never for production.
   run->SetStoreTraj(kFALSE);

   run->Init();

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);

   timer.Stop();
   std::cout << std::endl << std::endl;
   std::cout << "Macro finished succesfully." << std::endl;
   std::cout << "Output file is " << outFile << std::endl;
   std::cout << "Parameter file is " << parFile << std::endl;
   std::cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << "s" << std::endl << std::endl;
}
