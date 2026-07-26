/// @file fitUKF_C14.C
/// @brief Run ONLY the new UKF fitter on a pre-reconstructed a1954 14C(p,p') file.
///
/// Ported from a1975 UKF/pipeline/fitUKF_a1975.C. Reads <fileName>_reco.root
/// (AtPatternEvent, from unpackReco_C14.C) and runs EventFit::AtFitterUKF wrapped
/// in AtFitterUKFMulti. Fast to iterate (skips unpack/PSA/PRA).
///
/// HANDEDNESS: experimental data has the opposite helix handedness to simulation,
/// so Bz is flipped for exp data. bFieldSign = -1 (EXPERIMENTAL, default), +1 = sim.
/// If fits diverge / give unphysical KE, try the other sign first (validated on a1975).
///
/// gasDensity default 6.5e-5 g/cm^3 = H2 at 600 torr, ~293 K (a1954 target).
/// (a1975 used 9.0e-5 for H2 at 1 bar; scaled by 600/760.)  << CONFIRM against run conditions.
///
///   root -b -q 'fitUKF_C14.C("run_0055", -1, "proton", -1)'

#include <map>

namespace {
struct PSpec {
   double massMeV;
   int Z;
   int A;
};
const std::map<TString, PSpec> kParticleTable = {
   {"proton", {938.27208816, 1, 1}}, {"deuteron", {1875.61294257, 1, 2}}, {"triton", {2808.92113298, 1, 3}},
   {"He3", {2808.39160743, 2, 3}},   {"alpha", {3727.3794066, 2, 4}},
};
const double kU_MeV = 931.49410242;
const double kE_C = 1.602176634e-19;

std::unique_ptr<EventFit::AtFitterUKF> MakeHypothesis(const TString &name, int bFieldSign, double bFieldMag,
                                                      double gasDensity, double measSigma, double momSigmaFrac,
                                                      int nIter, int minClusters, bool refTrack = false,
                                                      double refInflation = 4.0)
{
   auto it = kParticleTable.find(name);
   if (it == kParticleTable.end()) {
      std::cout << "\033[1;31m[fitUKF_C14] Unknown particle '" << name << "' — skipping.\033[0m" << std::endl;
      return nullptr;
   }
   const PSpec &p = it->second;
   const double massAmu = p.massMeV / kU_MeV;
   const double charge = p.Z * kE_C;

   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(p.A, p.Z, massAmu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen target gas
   eloss->SetMaterial(mat);

   auto ukf = std::make_unique<EventFit::AtFitterUKF>(charge, p.massMeV, std::move(eloss));
   ukf->SetBField(ROOT::Math::XYZVector(0, 0, bFieldSign * bFieldMag));
   ukf->SetUKFParameters(1e-3, 2.0, 0.0);
   ukf->SetMeasurementSigma(measSigma);
   ukf->SetMomentumSigmaFrac(momSigmaFrac);
   ukf->SetEnableEnergyStraggling(false);
   ukf->SetMinClusters(minClusters);
   ukf->SetNIterations(nIter);
   ukf->SetUseRefTrack(refTrack);
   ukf->SetRefTrackInflation(refInflation);
   ukf->SetZPadPlane(1000.0);
   return ukf;
}
} // namespace

void fitUKF_C14(TString fileName = "run_0055", Long64_t nEvents = -1, TString particles = "proton",
                 Int_t bFieldSign = -1, Double_t bFieldMag = 2.85, Double_t gasDensity = 6.5e-5, TString outSuffix = "",
                 TString ioDir = "/home/yassid/a1954_C14_reco/", Double_t measSigma = 0.5, Double_t momSigmaFrac = 0.1,
                 Int_t nIter = 1, Int_t minClusters = 10, TString outDir = "", Bool_t refTrack = false,
                 Double_t refInflation = 4.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_ukf" + outSuffix + ".root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954_C14.par";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found. Run unpackReco_C14.C first.\033[0m" << std::endl;
      return;
   }

   std::cout << "\033[1;33m=== fitUKF_C14 ===\033[0m" << std::endl;
   std::cout << "  input      : " << inputFile << std::endl;
   std::cout << "  particles  : " << particles << std::endl;
   std::cout << "  bFieldSign : " << bFieldSign << "  (" << (bFieldSign < 0 ? "EXPERIMENTAL" : "simulation")
             << " handedness),  Bz = " << bFieldSign * bFieldMag << " T" << std::endl;
   std::cout << "  gasDensity : " << gasDensity << " g/cm^3" << std::endl;

   TStopwatch timer;

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
      auto hypo = MakeHypothesis(pname, bFieldSign, bFieldMag, gasDensity, measSigma, momSigmaFrac, nIter, minClusters,
                                 refTrack, refInflation);
      if (hypo) {
         multi->AddHypothesis(std::move(hypo), pname.Data(), 0);
         std::cout << "  + hypothesis: " << pname << std::endl;
      }
   }
   delete toks;
   if (multi->GetNHypotheses() == 0) {
      std::cout << "\033[1;31mNo valid particle hypotheses — aborting.\033[0m" << std::endl;
      return;
   }

   AtFitterTask *fitterTask = new AtFitterTask(std::move(multi));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   fRun->AddTask(fitterTask);
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents < 0 ? 0 : nEvents);
   timer.Stop();

   std::cout << "\n\033[1;32mDone.\033[0m  Output: " << outputFile << std::endl;
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}
