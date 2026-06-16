/// @file fitUKF_a1975_deuterium.C
/// @brief Run the OpenKF UKF fitter on the a1975 D2-target 16C(d,p)17C recos, PROTON
/// hypothesis, for a head-to-head comparison vs the genfit production. Mirrors
/// UKF/pipeline/fitUKF_a1975.C but: (1) reads reco_d2/ + the deuterium par, (2) adds
/// an AtPIDTask so the output carries AtPIDEvent (the SAME Spyral plane the proton
/// gate proton_band_d2_v2.json was built on), and (3) writes
/// <run>_genfitter_p_UKF.root so ex_dp_a1975.C can read it directly via fitSuffix="_UKF".
///
///   root -b -q 'fitUKF_a1975_deuterium.C("run_0016", -1, "proton", -1, 2.85, 9.0e-5, "/mnt/f/a1975/reco_d2/")'
///
/// bFieldSign=-1 = experimental handedness (default), same as the genfit B=-2.85.

#include <map>

namespace {
struct PSpecU {
   double massMeV;
   int Z;
   int A;
};
const std::map<TString, PSpecU> kParticleTableU = {
   {"proton", {938.27208816, 1, 1}}, {"deuteron", {1875.61294257, 1, 2}}, {"triton", {2808.92113298, 1, 3}}};
const double kU_MeVu = 931.49410242;
const double kE_Cu = 1.602176634e-19;

std::unique_ptr<EventFit::AtFitterUKF> MakeHypoU(const TString &name, int bFieldSign, double bFieldMag,
                                                 double gasDensity, double measSigma, double momSigmaFrac, int nIter,
                                                 int minClusters)
{
   auto it = kParticleTableU.find(name);
   if (it == kParticleTableU.end())
      return nullptr;
   const PSpecU &p = it->second;
   const double massAmu = p.massMeV / kU_MeVu;
   const double charge = p.Z * kE_Cu;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(p.A, p.Z, massAmu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1));
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
   return ukf;
}
} // namespace

void fitUKF_a1975_deuterium(TString fileName = "run_0016", Long64_t nEvents = -1, TString particles = "proton",
                            Int_t bFieldSign = -1, Double_t bFieldMag = 2.85, Double_t gasDensity = 9.0e-5,
                            TString ioDir = "/mnt/f/a1975/reco_d2/", TString outDir = "", Double_t measSigma = 2.0,
                            Double_t momSigmaFrac = 0.3, Int_t nIter = 1, Int_t minClusters = 10)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter_p_UKF.root"; // ex_dp reads via fitSuffix="_UKF"
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitUKF_a1975_deuterium (16C(d,p)17C, UKF, PROTON hyp) ===\033[0m\n"
             << "  in : " << inputFile << "\n  out: " << outputFile << "\n  Bz=" << bFieldSign * bFieldMag
             << " T  particles=" << particles << "\n";

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
      auto hypo = MakeHypoU(pname, bFieldSign, bFieldMag, gasDensity, measSigma, momSigmaFrac, nIter, minClusters);
      if (hypo) {
         multi->AddHypothesis(std::move(hypo), pname.Data(), 0);
         std::cout << "  + hypothesis: " << pname << "\n";
      }
   }
   delete toks;
   if (multi->GetNHypotheses() == 0) {
      std::cout << "\033[1;31mNo valid hypotheses.\033[0m\n";
      return;
   }

   // PID task so the output carries AtPIDEvent (gate plane), like the genfit macro.
   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(multi));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(pidTask);
   fRun->AddTask(fitterTask);

   TStopwatch timer;
   timer.Start();
   fRun->Init();
   fRun->Run(0, nEvents < 0 ? 0 : nEvents);
   timer.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << timer.RealTime() << " s)\n";
}
