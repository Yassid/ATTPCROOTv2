/// @file unpackNFit_a1975_UKF.C
/// @brief Full a1975 16C+p pipeline with the NEW UKF, end-to-end from raw HDF5.
///
///   raw HDF5 -> AtRawEvent -> PSAMax(AtEventH) -> SC(AtEventCorrected)
///            -> PRA(AtPatternEvent) -> UKF(AtTrackingEvent + AtFitMetadata)
///
/// This is the production entry point. For fast fitter iteration on an already
/// reconstructed file, use fitUKF_a1975.C (reads <run>_reco.root) instead.
///
/// PARTICLE SELECTION and ORIENTATION/HANDEDNESS work exactly as in
/// fitUKF_a1975.C — see the long comment there. Key point: experimental a1975
/// data has the OPPOSITE helix handedness to simulation (sim reverses z in
/// digitization, experiment does not), so the UKF uses Bz < 0 for exp data.
/// Empirically confirmed on run_0116: Bz=-2.85 gives chi2/ndf~0.04 and physical
/// proton KE, while Bz=+2.85 gives chi2/ndf~3 and unphysical KE.
///
/// Run: root -b -q 'unpackNFit_a1975_UKF.C("run_0116", 1000, "proton", -1)'

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
                                                      double gasDensity)
{
   auto it = kParticleTable.find(name);
   if (it == kParticleTable.end()) {
      std::cout << "\033[1;31m[UKF] Unknown particle '" << name << "' — skipping.\033[0m" << std::endl;
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
   ukf->SetBField(ROOT::Math::XYZVector(0, 0, bFieldSign * bFieldMag)); // sign = handedness
   ukf->SetUKFParameters(1e-3, 2.0, 0.0);
   ukf->SetMeasurementSigma(2.0);
   ukf->SetMomentumSigmaFrac(0.3);
   ukf->SetEnableEnergyStraggling(false);
   ukf->SetMinClusters(10);
   ukf->SetNIterations(1);
   ukf->SetZPadPlane(1000.0);
   return ukf;
}
} // namespace

void unpackNFit_a1975_UKF(TString fileName = "run_0116", Long64_t nEvents = 1000, TString particles = "proton",
                          Int_t bFieldSign = -1, Double_t bFieldMag = 2.85, Double_t gasDensity = 9.0e-5)
{
   gSystem->Load("libAtReconstruction.so");

   TStopwatch timer;
   timer.Start();

   TString parameterFile = "ATTPC.a1954.par";
   TString filepath = "/mnt/f/a1975/h5/";
   TString inputFile = filepath + fileName + ".h5";
   TString scriptfile = "ANL2023.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString outputFile = fileName + "_ukf_full.root";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString zlutFile = dir + "/resources/corrections/a1954/zLUT.txt";
   TString radlutFile = dir + "/resources/corrections/a1954/radLUT.txt";
   TString tralutFile = dir + "/resources/corrections/a1954/traLUT.txt";

   std::cout << "\033[1;33mInput : " << inputFile << "\033[0m" << std::endl;
   std::cout << "  particles=" << particles << "  bFieldSign=" << bFieldSign << " (Bz=" << bFieldSign * bFieldMag
             << " T)" << std::endl;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: input HDF5 not found.\033[0m" << std::endl;
      return;
   }

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   // --- Unpack ---
   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(false);

   // --- PSA ---
   auto psa = new AtPSAMax();
   psa->SetThreshold(20);
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(false);
   psaTask->SetOutputBranch("AtEventH");

   // --- Space-charge correction ---
   auto SCModel = std::make_unique<AtEDistortionModel>();
   SCModel->SetCorrectionMaps(zlutFile.Data(), radlutFile.Data(), tralutFile.Data());
   auto SCTask = new AtSpaceChargeCorrectionTask(std::move(SCModel));
   SCTask->SetInputBranch("AtEventH");

   // --- PRA via dependency injection (AtPRAtask takes the algorithm as a unique_ptr) ---
   // AtTrackFinderTC carries the standard TC defaults; only cluster radius/distance set.
   auto praAlgo = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   praAlgo->SetClusterRadius(15.0);
   praAlgo->SetClusterDistance(7.5);
   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch("AtEventCorrected");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);

   // --- UKF (multi-hypothesis from the particle list) ---
   auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
   TObjArray *toks = particles.Tokenize(",");
   for (int i = 0; i < toks->GetEntries(); ++i) {
      TString pname = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
      auto hypo = MakeHypothesis(pname, bFieldSign, bFieldMag, gasDensity);
      if (hypo) {
         multi->AddHypothesis(std::move(hypo), pname.Data(), 0);
         std::cout << "  + UKF hypothesis: " << pname << std::endl;
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
   fitterTask->SetPersistence(true);

   run->AddTask(unpackTask);
   run->AddTask(psaTask);
   run->AddTask(SCTask);
   run->AddTask(praTask);
   run->AddTask(fitterTask);

   std::cout << "***** Init ******" << std::endl;
   run->Init();

   auto numEvents = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < numEvents)
      numEvents = nEvents;
   std::cout << "\033[1;32mProcessing " << numEvents << " events end-to-end (unpack..UKF).\033[0m" << std::endl;

   run->Run(0, numEvents);

   std::cout << "\nDone. Output: " << outputFile << std::endl;
   timer.Stop();
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}
