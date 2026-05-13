/// @file run_ukf_HYDRA.C
/// @brief Run the UKF fitter on digitized HYDRA pion events.
///
/// Reads ./data/output_digi<inSuffix>.root and writes
/// ./data/output_ukf_HYDRA<outSuffix>.root. Mirrors pi_TPC's run_ukf_only.C
/// with the HYDRA geometry (drift length 150 mm along z, beam axis along
/// ATTPCROOT-y so the back-extrapolation surface is y=-50 cm).
///
/// Run: root -b -q 'run_ukf_HYDRA.C(1000, +1)'

void run_ukf_HYDRA(int nEvents = 10000, int pionSign = +1,
                   double momSigmaFrac = 0.1, int nIterations = 3,
                   const char *outSuffix = "", double eLossScale = 1.0,
                   const char *inSuffix = "", double Bz_T = 2.0)
{
   FairLogger::GetLogger()->SetLogScreenLevel("INFO");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi" + inSuffix + ".root";
   TString outputFile = inOutDir + "output_ukf_HYDRA" + outSuffix + ".root";
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
   // HYDRA active z (drift) spans 0–150 mm in ATTPCROOT-z. Pad plane is at
   // z=0 in lab; the digi convention puts it at SetZPadPlane (set so
   // z_lab = ZPadPlane − Z_digi). Use the drift length (150 mm) for the
   // back-extrap cap since the vertex is at z_lab ≈ 75 mm (mid-drift).
   ukfFitter->SetZPadPlane(150.0);
   ukfFitter->SetBackExtrapMaxPath(250.0);
   ukfFitter->SetUseClusterDirSeed(true);
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
   std::cout << "B = " << Bz_T << " T  /  KE seed range from PRA Brho\n";
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
