/// @file run_ukf_multi.C
/// @brief Multi-particle-hypothesis UKF fit for PUMA events.
///        Each PRA track is fitted under K+, K-, pi+, pi-; the best
///        reduced chi^2/ndf wins. B = 4 T, gas = P10 @ 1 bar.
///
/// Run: root -b -q 'run_ukf_multi.C(50)'

void run_ukf_multi(int nEvents = 50, TString tag = "")
{
   FairLogger::GetLogger()->SetLogScreenLevel("INFO");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_digi" + tag + ".root";
   TString outputFile = inOutDir + "output_ukf_multi" + tag + ".root";
   TString paramFile = "ATTPC.PUMA_sim.par";

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

   const double e_C = 1.602176634e-19;
   const double u_MeV = 931.49410372;
   const double m_K_MeV = 493.677;
   const double m_pi_MeV = 139.57039;

   const double masses[4] = {m_K_MeV, m_K_MeV, m_pi_MeV, m_pi_MeV};
   const int signs[4] = {+1, -1, +1, -1};
   const std::string names[4] = {"K+", "K-", "pi+", "pi-"};

   auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();

   for (int i = 0; i < 4; ++i) {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
      eloss->SetProjectile(1, 1, masses[i] / u_MeV);
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back(std::make_tuple(18, 40, 9));
      mat.push_back(std::make_tuple(6, 12, 1));
      mat.push_back(std::make_tuple(1, 1, 4));
      eloss->SetMaterial(mat);

      auto ukf = std::make_unique<EventFit::AtFitterUKF>(signs[i] * e_C, masses[i], std::move(eloss));
      ukf->SetBField(ROOT::Math::XYZVector(0, 0, 4.0));
      ukf->SetUKFParameters(1e-3, 2.0, 0.0);
      ukf->SetMeasurementSigma(1.0);
      ukf->SetMomentumSigmaFrac(0.1);
      ukf->SetEnableEnergyStraggling(true);
      ukf->SetMinClusters(2);
      ukf->SetNIterations(3);
      ukf->SetZPadPlane(0.0);
      ukf->SetBackExtrapMaxPath(250.0);
      multi->AddHypothesis(std::move(ukf), names[i]);
   }

   AtFitterTask *fitterTask = new AtFitterTask(std::move(multi));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(fitterTask);
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "\nMulti-hypothesis UKF done.  Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n";
}
