/// @file run_ukf_only.C
/// @brief Run the UKF fitter on digitized pion events.
///
/// Reads ./data/output_digi.root and writes ./data/output_ukf_only.root.
/// Mass / charge / dE/dx are configured for charged pions.
///
/// CATIMA caveat: AtELossCATIMA is built around nuclear projectiles
/// parameterised by (A, Z) plus a mass-in-amu. We pass A=1, Z=1 (so the
/// Bethe-Bloch Z^2/beta^2 scaling is correct for a unit-charge particle)
/// and mass_amu = m_pi/u so the kinetic-energy-per-amu argument used
/// internally yields the right beta. This is a first-order approximation:
/// shell corrections, density effect and nuclear inelastic processes are
/// not pion-specific. Adequate for UKF kinematics tests; revisit if you
/// need few-percent dE/dx accuracy.
///
/// Run: root -b -q 'run_ukf_only.C(10000, +1)'   for pi+
///      root -b -q 'run_ukf_only.C(10000, -1)'   for pi-
///      root -b -q 'run_ukf_only.C(2000, +1, 0.3, 3, "_mom03")'  override

void run_ukf_only(int nEvents = 10000, int pionSign = +1,
                  double momSigmaFrac = 0.1, int nIterations = 3,
                  const char *outSuffix = "", double eLossScale = 1.0,
                  const char *inSuffix = "", double Bz_T = 2.0)
{
   // INFO captures the per-track UKF rejection reasons (the four
   // "AtFitterUKF: ... Skipping." messages, see below).
   FairLogger::GetLogger()->SetLogScreenLevel("INFO");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi" + inSuffix + ".root";
   TString outputFile = inOutDir + "output_ukf_only" + outSuffix + ".root";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // ---- pion-specific parameters --------------------------------------
   const double e_C = 1.602176634e-19;       // Coulombs
   const double mass_pi_MeV = 139.57039;     // PDG charged pion mass
   const double u_MeV = 931.49410372;        // 1 amu in MeV/c^2
   const double mass_pi_amu = mass_pi_MeV / u_MeV;

   const double charge = (pionSign >= 0) ? +e_C : -e_C;

   // P10 (90% Ar + 10% CH4) at 1 bar, 273 K -> rho = 1.654e-3 g/cm^3
   //
   // TODO: replace with a proper pion energy-loss model. CATIMA is a nuclear-
   // projectile model; passing A=Z=1 with mass_amu = m_pi/u recovers Bethe-Bloch
   // Z^2/beta^2 scaling but misses pion-specific shell/density corrections and
   // nuclear inelastic. Adequate for relativistic MIPs (few-percent dE/dx);
   // revisit before trusting Bragg-curve / low-energy stopping-pion work.
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
   eloss->SetProjectile(1, 1, mass_pi_amu); // CATIMA workaround, see header
   // Material composition for CATIMA: (Z, A, n_atoms_per_molecule).
   // Use the elemental decomposition: 0.9 Ar (atom-fraction) + 0.1 CH4
   // CATIMA treats this as a fixed compound; 1 Ar + 0.1 CH4 ~ approximated as
   // (Ar) + (C) + (4*0.1/0.9 H) keeping Ar-dominated stopping behavior.
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(18, 40, 9));  // Ar  (9 in P10 per "unit")
   mat.push_back(std::make_tuple(6, 12, 1));   // C   (1 from CH4)
   mat.push_back(std::make_tuple(1, 1, 4));    // H   (4 from CH4)
   eloss->SetMaterial(mat);

   auto ukfFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass_pi_MeV, std::move(eloss));
   ukfFitter->SetBField(ROOT::Math::XYZVector(0, 0, Bz_T)); // match pi_TPC_sim.C
   ukfFitter->SetUKFParameters(1e-3, 2.0, 0.0);
   ukfFitter->SetMeasurementSigma(1.0);
   ukfFitter->SetMomentumSigmaFrac(momSigmaFrac);
   ukfFitter->SetEnableEnergyStraggling(true);
   ukfFitter->SetELossScaleFactor(eLossScale);
   ukfFitter->SetMinClusters(5);
   ukfFitter->SetNIterations(nIterations);
   ukfFitter->SetZPadPlane(1000.0);
   // HELIOS-style: vertex deep inside the TPC, first cluster can be ~80-150 mm
   // away — open up the back-extrapolation cap (default 50 mm).
   ukfFitter->SetBackExtrapMaxPath(250.0);
   // Circle-tangent seed (vs legacy chord seed). The chord seed is wrong for
   // backward (theta_lab>90°) tracks because seed=highest-Z_digi cluster is
   // the track *end*, so the chord points opposite to motion. The tangent
   // logic sorts clusters by xy-distance from the beam axis and uses the
   // PRA circle tangent — direction-agnostic and works isotropically.
   ukfFitter->SetUseClusterDirSeed(true);
   // Arc-walk re-clustering: geometry-based, gap-immune, targets a fixed
   // number of well-placed clusters along the helix arc. Better than
   // Smooth3D's distance-based clustering for long MIP tracks where pad
   // gaps and diffusion can leave Smooth3D with under/oversized clusters.
   ukfFitter->SetUseArcWalk(true);
   ukfFitter->SetTargetClusters(15);
   ukfFitter->SetAdaptiveDistBounds(4.0, 14.0);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(ukfFitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(fitterTask);
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\nMacro finished succesfully.\n";
   std::cout << "PDG: " << ((pionSign >= 0) ? "+211 (pi+)" : "-211 (pi-)") << "\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
