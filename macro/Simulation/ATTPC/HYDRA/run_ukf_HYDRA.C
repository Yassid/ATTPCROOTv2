/// @file run_ukf_HYDRA.C
/// @brief Run the UKF fitter on digitized HYDRA Prototype pion events.
///
/// Reads ./data/output_digi<inSuffix>.root and writes
/// ./data/output_ukf_HYDRA<outSuffix>.root. Mirrors pi_TPC's run_ukf_only.C
/// with the HYDRA Prototype geometry (drift length 294 mm along +z;
/// beam axis is ATTPCROOT-x = HYDRA-z; lower-left of active region at
/// world (0, 0, 0)).
///
/// Run: root -b -q 'run_ukf_HYDRA.C(1000, -1)'

void run_ukf_HYDRA(int nEvents = 10000, int pionSign = -1,
                   double momSigmaFrac = 0.1, int nIterations = 3,
                   const char *outSuffix = "", double eLossScale = 1.0,
                   const char *inSuffix = "", double Bz_T = 2.0)
{
   FairLogger::GetLogger()->SetLogScreenLevel("INFO");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi" + inSuffix + ".root";
   TString outputFile = inOutDir + "output_ukf_HYDRA" + outSuffix + ".root";
   TString paramFile = "ATTPC.HYDRA_sim.par";

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
   const double e_C = 1.602176634e-19;
   const double mass_pi_MeV = 139.57039;
   const double u_MeV = 931.49410372;
   const double mass_pi_amu = mass_pi_MeV / u_MeV;
   const double charge = (pionSign >= 0) ? +e_C : -e_C;

   // P10 at 1 bar — same gas as pi_TPC.
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
   eloss->SetProjectile(1, 1, mass_pi_amu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(18, 40, 9));
   mat.push_back(std::make_tuple(6, 12, 1));
   mat.push_back(std::make_tuple(1, 1, 4));
   eloss->SetMaterial(mat);

   auto ukfFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass_pi_MeV, std::move(eloss));
   ukfFitter->SetBField(ROOT::Math::XYZVector(0, 0, Bz_T));
   ukfFitter->SetUKFParameters(1e-3, 2.0, 0.0);
   ukfFitter->SetMeasurementSigma(1.0);
   ukfFitter->SetMomentumSigmaFrac(momSigmaFrac);
   ukfFitter->SetEnableEnergyStraggling(true);
   ukfFitter->SetELossScaleFactor(eLossScale);
   ukfFitter->SetMinClusters(5);
   ukfFitter->SetNIterations(nIterations);
   // HYDRA Prototype active z (drift) spans 0–294 mm in ATTPCROOT-z. Pad
   // plane is at z=0 in lab; the digi convention puts it at SetZPadPlane
   // (z_lab = ZPadPlane − Z_digi). Vertex at (-40, 4.4, 14.7) cm =
   // z_lab ≈ 147 mm (mid-drift). Back-extrap cap covers vertex → entrance
   // (~28 cm in beam direction) plus margin.
   ukfFitter->SetZPadPlane(294.0);
   ukfFitter->SetBackExtrapMaxPath(450.0);
   ukfFitter->SetUseClusterDirSeed(true);
   ukfFitter->SetUseArcWalk(true);
   // Use the helix-POCA back-extrap (closed-form circle + helix-pitch z)
   // instead of the linear-tangent default. Required for HYDRA's visibly
   // curved low-p tracks (sagitta ~25 mm at 200 MeV/c over 256 mm chord
   // is too large for the straight-line approximation) AND so the
   // φ-rotation step below can actually fire.
   ukfFitter->SetUseHelixBackExtrap(true);
   // Make Kinematics.phi the angle at the back-extrapolated vertex, not
   // at the first cluster. Removes the first-cluster anchor wobble that
   // inflates σ_φ at low p (38 mrad → ~13 mrad at 200 MeV/c).
   ukfFitter->SetUpdateAnglesOnBackExtrap(true);
   // Straight-line tail disabled: the widened field region in HYDRA_sim
   // now covers the target plane upstream, so the pion curves along a
   // helix from production vertex all the way to the chamber — there's
   // no field-free region to straight-extrap through. The helix POCA
   // endpoint of the UKF (~chamber face) is the right place to stop.
   // ukfFitter->SetBackExtrapTargetX(-400.0);
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
   std::cout << "B = " << Bz_T << " T  /  KE seed range from PRA Brho\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
