/// @file fitBoth_a1975.C
/// @brief Fit a pre-reconstructed a1975 16C+p file with BOTH Kalman filters in one
/// pass, writing two separate tracking branches so the event display can overlay them:
///     AtTrackingEventGenfit  — EventFit::AtGenfitter (genfit KalmanFitterRefTrack)
///     AtTrackingEventUKF      — EventFit::AtFitterUKF  (unscented KF + CATIMA)
/// The source hits (AtEventCorrected) and AtPatternEvent are carried forward, so the
/// output file is self-contained for run_eve_both_a1975.C.
///
/// Both fitters use the SAME proton hypothesis (the validated 16C(p,p) channel) so the
/// two polylines are an apples-to-apples comparison. Keep nEvents modest for display.
///
///   root -b -q 'pipeline/fitBoth_a1975.C("run_0106", 200, "/mnt/f/a1975/reco/", "/tmp/")'

void fitBoth_a1975(TString fileName = "run_0106", Long64_t nEvents = 200, TString ioDir = "/mnt/f/a1975/reco/",
                   TString outDir = "", Int_t bFieldSign = -1, Double_t bFieldMag = 2.85,
                   TString pidGate = "pid/proton_band.json", Double_t measSigma = 4.0,
                   Double_t gasDensity = 9.0e-5)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   const double kU_MeV = 931.49410242, kE_C = 1.602176634e-19;
   // proton hypothesis (shared by both fitters)
   const Int_t pdg = 2212;
   const Double_t protonMassMeV = 938.27208816;
   const Double_t massAmu = protonMassMeV / kU_MeV; // 1.00727...
   const Int_t Z = 1, A = 1;

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_both.root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitBoth_a1975 (genfit + UKF) ===\033[0m\n  in  : " << inputFile
             << "\n  out : " << outputFile << "\n  Bz = " << bFieldSign * bFieldMag << " T\n";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   if (gROOT->FindObject("FAIRGeom") == nullptr) { // genfit needs gGeoManager
      TFile *gf = TFile::Open(geoManFile);
      gf->Get("FAIRGeom");
   }

   // ─── genfit (signed B for experimental handedness; internal Bethe-Bloch, mat off) ───
   auto genfitter = std::make_unique<EventFit::AtGenfitter>(bFieldSign * bFieldMag, pdg, massAmu, Z,
                                                            /*eLossFile*/ "", /*noMatEffects*/ kTRUE, 2, 5);
   genfitter->SetZPadPlane(1000.0);
   genfitter->SetMeasSigma(measSigma);
   if (pidGate.Length() && !gSystem->AccessPathName(pidGate.Data()))
      genfitter->SetPIDGate(pidGate.Data());
   auto genfitTask = new AtFitterTask(std::move(genfitter));
   genfitTask->SetInputBranch("AtPatternEvent");
   genfitTask->SetOutputBranch("AtTrackingEventGenfit");
   genfitTask->SetFitMetadataBranch("AtFitMetadataGenfit");
   genfitTask->SetPersistence(kTRUE);

   // ─── UKF (proton hypothesis; same calibration as fitUKF_a1975.C) ───
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(A, Z, massAmu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen target gas
   eloss->SetMaterial(mat);
   auto ukf = std::make_unique<EventFit::AtFitterUKF>(Z * kE_C, protonMassMeV, std::move(eloss));
   ukf->SetBField(ROOT::Math::XYZVector(0, 0, bFieldSign * bFieldMag));
   ukf->SetUKFParameters(1e-3, 2.0, 0.0);
   ukf->SetMeasurementSigma(2.0);
   ukf->SetMomentumSigmaFrac(0.3);
   ukf->SetEnableEnergyStraggling(false);
   ukf->SetMinClusters(5); // 5 (not 10) recovers short tracks the viewer showed UKF dropping
                           // (+29% yield 10->4, physical-KE frac 71->68%; saturates at 4 — ukfMinClus scan)
   ukf->SetNIterations(1);
   ukf->SetZPadPlane(1000.0);
   auto ukfTask = new AtFitterTask(std::move(ukf));
   ukfTask->SetInputBranch("AtPatternEvent");
   ukfTask->SetOutputBranch("AtTrackingEventUKF");
   ukfTask->SetFitMetadataBranch("AtFitMetadataUKF");
   ukfTask->SetPersistence(kTRUE);

   run->AddTask(genfitTask);
   run->AddTask(ukfTask);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   t.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << t.RealTime() << " s)\n"
             << "  branches: AtTrackingEventGenfit, AtTrackingEventUKF (+ AtEventCorrected, AtPatternEvent)\n";
}
