/// @file fitUKF_C15d.C
/// @brief 15C + d  --  stage 2 alternative: OpenKF UKF fit of a C15d reco.
///
/// Same input and the same PID plane as fitGenfit_C15d.C, so the two are a controlled
/// comparison on ONE variable (the fitter) over identical reconstruction. Writes
/// <run>_ukf_<tag><suffix>.root with AtTrackingEvent + AtPIDEvent.
///
///   root -b -q 'fitUKF_C15d.C("run_0017", 500, "proton", "/home/yassid/C15d_reco/")'
///
/// The UKF gets its energy loss from CATIMA directly (AtTools::AtELossCATIMA), so with
/// genfit's CATIMA backend on, both fitters describe the gas with the SAME physics and the
/// same density -- which is the only way the comparison isolates the estimator.
///
/// matA = 2 and the (2,1,1) target tuple: the target is DEUTERIUM. The default elsewhere is
/// ordinary hydrogen (1,1,1), which has the same Z but half the A -- twice the electrons per
/// unit mass, so twice the stopping per g/cm2. Getting this wrong is a ~2x error in dE/dx
/// that no fit-quality number will flag.
///
/// BACKWARD TRACKS: in (d,p) inverse kinematics the proton goes largely backward, so
///   clusterDirSeed = true  -> seed direction from the cluster geometry, not the
///                             half-sphere-ambiguous PRA GeoTheta
///   minClusters    = 4     -> backward tracks are cluster-poor; the usual 10 drops most
/// Both are existing AtFitterUKF knobs; nothing in the framework is modified.

#include <map>

namespace {
struct PSpecU {
   double massMeV;
   int Z;
   int A;
};
const std::map<TString, PSpecU> kParticleTableC15d = {
   {"proton", {938.27208816, 1, 1}}, {"deuteron", {1875.61294257, 1, 2}}, {"triton", {2808.92113298, 1, 3}}};
const double kU_MeVu_C15d = 931.49410242;
const double kE_C_C15d = 1.602176634e-19;

std::unique_ptr<EventFit::AtFitterUKF> MakeHypoC15d(const TString &name, int bFieldSign, double bFieldMag, int matA,
                                                    double gasDensity, double measSigma, double momSigmaFrac,
                                                    int nIter, int minClusters, bool clusterDirSeed)
{
   auto it = kParticleTableC15d.find(name);
   if (it == kParticleTableC15d.end())
      return nullptr;
   const PSpecU &p = it->second;
   const double massAmu = p.massMeV / kU_MeVu_C15d;
   const double charge = p.Z * kE_C_C15d;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(p.A, p.Z, massAmu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(matA, 1, 1)); // D2: A=2, Z=1, one atom per formula unit
   eloss->SetMaterial(mat);
   auto ukf = std::make_unique<EventFit::AtFitterUKF>(charge, p.massMeV, std::move(eloss));
   ukf->SetBField(ROOT::Math::XYZVector(0, 0, bFieldSign * bFieldMag));
   ukf->SetUKFParameters(1e-3, 2.0, 0.0);
   ukf->SetMeasurementSigma(measSigma);
   ukf->SetMomentumSigmaFrac(momSigmaFrac);
   ukf->SetEnableEnergyStraggling(false);
   ukf->SetMinClusters(minClusters);
   ukf->SetNIterations(nIter);
   ukf->SetZPadPlane(1000.0);
   ukf->SetUseClusterDirSeed(clusterDirSeed);
   return ukf;
}
} // namespace

void fitUKF_C15d(TString fileName = "run_0017", Long64_t nEvents = -1, TString particles = "proton",
                 TString ioDir = "/home/yassid/C15d_reco/", TString outDir = "", Int_t bFieldSign = -1,
                 Double_t bFieldMag = 2.85, Double_t gasDensity = 6.5643e-5, Double_t measSigma = 2.0,
                 Double_t momSigmaFrac = 0.3, Int_t nIter = 1, Int_t minClusters = 4,
                 Bool_t clusterDirSeed = kTRUE, TString recoSuffix = "_reco", TString outSuffix = "",
                 TString parName = "ATTPC.C15d_a2091_D2.par", Int_t matA = 2, Bool_t gainMatch = kTRUE,
                 TString gainTable = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   if (dir.Length() == 0) {
      std::cout << "\033[1;31mERROR: VMCWORKDIR unset -- source build/config.sh first.\033[0m\n";
      return;
   }
   TString inputFile = ioDir + fileName + recoSuffix + ".root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString tag = particles;
   tag.ReplaceAll(",", "_");
   TString outputFile = outBase + fileName + "_ukf_" + tag + outSuffix + ".root";
   TString digiParFile = dir + "/parameters/" + parName;

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitUKF_C15d (15C + d, D2 300 torr, UKF, hyp: " << particles << ") ===\033[0m\n"
             << "  in  : " << inputFile << "\n  out : " << outputFile << "\n  par : " << parName
             << "\n  Bz=" << bFieldSign * bFieldMag << " T   rho=" << gasDensity << " g/cm3   target A=" << matA
             << "\n";

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);
   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
   TObjArray *toks = particles.Tokenize(",");
   for (int i = 0; i < toks->GetEntries(); ++i) {
      TString pname = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
      auto hypo = MakeHypoC15d(pname, bFieldSign, bFieldMag, matA, gasDensity, measSigma, momSigmaFrac, nIter,
                               minClusters, clusterDirSeed);
      if (hypo) {
         multi->AddHypothesis(std::move(hypo), pname.Data(), 0);
         std::cout << "  + hypothesis: " << pname << "\n";
      } else {
         std::cout << "\033[1;31m  unknown species '" << pname << "' -- skipped\033[0m\n";
      }
   }
   delete toks;
   if (multi->GetNHypotheses() == 0) {
      std::cout << "\033[1;31mNo valid hypotheses.\033[0m\n";
      return;
   }

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(multi));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   // Opt-in per-run gain matching, added AFTER AtPIDTask so it rescales that task's output.
   // The UKF itself never looks at dE/dx, so unlike the genfit macro there is no in-fit gate
   // to fall out of step with the matched plane.
   AtGainMatchTask *gainTask = nullptr;
   if (gainMatch) {
      TString tbl = gainTable.Length() ? gainTable : (dir + "/macro/Unpack_HDF5/C15d/gainmatch_C15d.csv");
      const Int_t runNo = AtGainMatchTask::RunNumberFromName(fileName);
      if (runNo < 0) {
         std::cout << "\033[1;31mERROR: cannot parse a run number from '" << fileName
                   << "' -- gain matching would use the wrong factor. Aborting.\033[0m\n";
         return;
      }
      gainTask = new AtGainMatchTask(tbl.Data(), runNo);
      std::cout << "  \033[1;32mGain match: ON  (run " << runNo << ", table " << tbl << ")\033[0m\n";
   }

   fRun->AddTask(pidTask);
   if (gainTask)
      fRun->AddTask(gainTask);
   fRun->AddTask(fitterTask);

   TStopwatch timer;
   timer.Start();
   fRun->Init();
   fRun->Run(0, nEvents < 0 ? 0 : nEvents);
   timer.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << timer.RealTime() << " s)\n";
}
