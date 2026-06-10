/// @file fitUKF_a1975.C
/// @brief Run ONLY the new UKF fitter on a pre-reconstructed a1975 16C+p file.
///
/// Reads <fileName>_reco.root (must contain AtPatternEvent, produced by
/// unpackReco_a1975_UKF.C) and runs the UKF (EventFit::AtFitterUKF, wrapped in
/// AtFitterUKFMulti for multi-particle hypotheses). Fast to iterate because it
/// skips the unpack/PSA/SC/PRA stages.
///
/// ─── PARTICLE SELECTION (macro level) ────────────────────────────────────
/// `particles` is a comma-separated list, e.g. "proton" or "proton,deuteron".
/// Each name becomes one UKF hypothesis; the multi-fitter keeps the best
/// reduced-chi2/ndf per track and records the chosen particle in the output.
/// Supported: proton, deuteron, triton, He3, alpha.
///
/// ─── ORIENTATION / HANDEDNESS (critical for EXPERIMENTAL data) ────────────
/// In SIMULATION the digitization reverses z:  z_digi = ZPadPlane - z_MC
/// (AtClusterize.cxx:128). SetZPadPlane(1000) inverts this to recover the lab
/// frame, and the PRA charge-sign + UKF helix handedness are calibrated to that
/// SIM convention with B = +2.85 ẑ.
///
/// EXPERIMENTAL data never passes through that digitization reversal — its z
/// comes straight from drift time in PSA — so the spiral handedness is the
/// OPPOSITE of simulation. To make the UKF's propagated helix match the data we
/// flip the sign of Bz for experimental data. That is controlled by `bFieldSign`:
///     bFieldSign = -1  → experimental a1975 data (DEFAULT)
///     bFieldSign = +1  → simulation convention
/// If fits diverge / yield unphysical charges, try the other sign — this is the
/// first thing to validate. (A z-reflection ≡ Bz sign flip for helix handedness.)
///
/// Run: root -b -q 'fitUKF_a1975.C("run_0116", 1000, "proton", -1)'

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

/// Build one fully-configured UKF hypothesis for a named particle.
std::unique_ptr<EventFit::AtFitterUKF> MakeHypothesis(const TString &name, int bFieldSign, double bFieldMag,
                                                      double gasDensity, double measSigma, double momSigmaFrac,
                                                      int nIter, int minClusters)
{
   auto it = kParticleTable.find(name);
   if (it == kParticleTable.end()) {
      std::cout << "\033[1;31m[fitUKF_a1975] Unknown particle '" << name << "' — skipping.\033[0m" << std::endl;
      return nullptr;
   }
   const PSpec &p = it->second;
   const double massAmu = p.massMeV / kU_MeV;
   const double charge = p.Z * kE_C; // magnitude; handedness handled by Bz sign

   // Energy loss in the H 1 bar target gas (CATIMA; A,Z,mass of the projectile).
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(p.A, p.Z, massAmu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen target gas
   eloss->SetMaterial(mat);

   auto ukf = std::make_unique<EventFit::AtFitterUKF>(charge, p.massMeV, std::move(eloss));
   // a1975 16C+p legacy UKF path (matches macro/Simulation/ATTPC/16C_pp/run_reco_ukf.C),
   // with the Bz sign flipped for experimental-data handedness.
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

void fitUKF_a1975(TString fileName = "run_0116", Long64_t nEvents = -1, TString particles = "proton",
                  Int_t bFieldSign = -1, Double_t bFieldMag = 2.85, Double_t gasDensity = 9.0e-5,
                  TString outSuffix = "", TString ioDir = "", Double_t measSigma = 2.0, Double_t momSigmaFrac = 0.3,
                  Int_t nIter = 1, Int_t minClusters = 10, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = ioDir + fileName + "_reco.root";
   // outDir (if given) keeps fit outputs separate from the input/reco dir — e.g.
   // the (p,d) deuteron fits go to /mnt/f/a1975/reco_pd/ while reading reco from
   // /mnt/f/a1975/reco/. Falls back to ioDir when empty.
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_ukf" + outSuffix + ".root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found. Run unpackReco_a1975_UKF.C first.\033[0m"
                << std::endl;
      return;
   }

   std::cout << "\033[1;33m=== fitUKF_a1975 ===\033[0m" << std::endl;
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

   // Build the multi-hypothesis fitter from the requested particle list.
   auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
   TObjArray *toks = particles.Tokenize(",");
   for (int i = 0; i < toks->GetEntries(); ++i) {
      TString pname = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
      auto hypo = MakeHypothesis(pname, bFieldSign, bFieldMag, gasDensity, measSigma, momSigmaFrac, nIter, minClusters);
      if (hypo) {
         // chargeSign = 0: do NOT filter by PRA charge sign. The PRA sign is
         // calibrated on SIM handedness and is unreliable for experimental data
         // until the orientation is validated, so evaluate every hypothesis.
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
