/// @file UKFTask.C
/// @brief Demonstrates running EventFit::AtFitterUKF inside AtFitterTask.
///
/// Prerequisites:
///   - A digitized+PSA+PRA output file (AtPatternEvent branch) from a simulation
///     (e.g., produced by the ATTPC H2 simulation with proton tracking).
///   - VMCWORKDIR must point to the ATTPCROOTv2 source tree.
///   - The environment must be loaded: source build/config.sh
///
/// Run from the ROOT prompt:
///   root -l -q 'UKFTask.C("input_pra.root", "output_ukf.root")'

#include "AtELossCATIMA.h"
#include "AtFitterUKF.h"
#include "AtFitterTask.h"

#include <FairParAsciiFileIo.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

#include <TStopwatch.h>
#include <TString.h>

#include <memory>

void UKFTask(TString inputFile = "output_pra.root", TString outputFile = "output_ukf.root",
             TString paramFile = "")
{
   TStopwatch timer;
   timer.Start();

   TString dir = gSystem->Getenv("VMCWORKDIR");

   // ---- FairRunAna setup --------------------------------------------------
   FairRunAna *fRun = new FairRunAna();
   fRun->SetInputFile(inputFile);
   fRun->SetOutputFile(outputFile);

   if (paramFile.IsNull())
      paramFile = dir + "/parameters/AT.H2.par";

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(paramFile.Data(), "in");
   rtdb->setFirstInput(parIo);

   // ---- Build UKF fitter --------------------------------------------------
   // Proton (Z=1, mass=938.272 MeV/c^2) in 300 Torr H2 gas.
   constexpr double charge_p = 1.602176634e-19; // Coulombs
   constexpr double mass_p = 938.272;           // MeV/c^2

   // Energy-loss model: CATIMA for in-situ use (no file required at runtime).
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5); // density in g/cm^3
   eloss->SetProjectile(1, 1, 1);                                    // proton: A=1, Z=1, charge=1
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1}); // H2: A=1, Z=1, stoichiometry=1
   eloss->SetMaterial(mat);

   auto fitter = std::make_unique<EventFit::AtFitterUKF>(charge_p, mass_p, std::move(eloss));
   fitter->SetBField({0, 0, 2.85});       // 2.85 T solenoidal field along Z
   fitter->SetEField({0, 0, 0});          // no electric field correction
   fitter->SetUKFParameters(1e-3, 2, 0);  // alpha, beta, kappa
   fitter->SetMeasurementSigma(1.0);      // 1 mm position sigma
   fitter->SetMomentumSigmaFrac(0.1);     // 10% fractional momentum seed uncertainty
   fitter->SetMinClusters(3);
   fitter->SetEnableEnergyStraggling(true); // CATIMA supports GetRangeVariance()

   // ---- AtFitterTask ------------------------------------------------------
   AtFitterTask *fitterTask = new AtFitterTask(std::move(fitter));
   fitterTask->SetPersistence(kTRUE);
   // Branch names assume the default PRA output: "AtPatternEvent" → "AtTrackingEvent"
   // fitterTask->SetInputBranch("AtPatternEvent");   // default
   // fitterTask->SetOutputBranch("AtTrackingEvent"); // default

   fRun->AddTask(fitterTask);

   // ---- Init and run -------------------------------------------------------
   fRun->Init();
   fRun->Run(0, 0); // process all events

   timer.Stop();
   timer.Print();
}
