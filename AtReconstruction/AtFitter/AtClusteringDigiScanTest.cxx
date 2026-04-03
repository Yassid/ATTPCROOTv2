#include "AtELossCATIMA.h"
#include "AtFitMetadata.h"
#include "AtFittedTrack.h"
#include "AtFitterUKF.h"
#include "AtHitCluster.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"
#include "AtTrackingEvent.h"

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <FairLogger.h>

#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

// ===========================================================================
// Scan clustering parameters on real digitized data
// ===========================================================================
class ClusteringDigiScanTest : public testing::Test {
protected:
   static constexpr double kMass = 938.272;
   static constexpr double kCharge = 1.602176634e-19;
   static constexpr double kBz = 2.85;
   static constexpr double kGasDensity = 3.553e-5;

   struct FitResult {
      int event;
      double radius;
      double distance;
      int nClusters;
      double kineticEnergy;
      double theta;
      bool converged;
      bool good; // KE in range AND theta in range
   };

   std::unique_ptr<EventFit::AtFitterUKF> CreateFitter(bool enableStraggling = false)
   {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(kGasDensity);
      eloss->SetProjectile(1, 1, 1);
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back({1, 1, 1});
      eloss->SetMaterial(mat);

      auto fitter = std::make_unique<EventFit::AtFitterUKF>(kCharge, kMass, std::move(eloss));
      fitter->SetBField({0, 0, kBz});
      fitter->SetUKFParameters(1e-3, 2.0, 0.0);
      fitter->SetMeasurementSigma(2.0);
      fitter->SetMomentumSigmaFrac(0.3);
      fitter->SetEnableEnergyStraggling(enableStraggling);
      fitter->SetMinClusters(5);
      fitter->SetZPadPlane(1000.0);
      return fitter;
   }

   FitResult FitTrack(AtTrack &track, int eventId, bool straggling = false)
   {
      FitResult result{eventId, 0, 0, 0, -1, -1, false, false};

      result.nClusters = track.GetHitClusterArray()->size();
      if (result.nClusters < 5)
         return result;

      auto fitter = CreateFitter(straggling);
      fitter->Init();

      AtTrackingEvent trackingEvent;
      AtPatternEvent patternEvent;
      patternEvent.GetTrackCand().push_back(track);

      fitter->FitEvent(&trackingEvent, &patternEvent);

      auto &fittedTracks = trackingEvent.GetFittedTracks();
      if (fittedTracks.empty())
         return result;

      result.converged = true;
      auto kin = fittedTracks[0]->GetKinematics();
      result.kineticEnergy = kin.kineticEnergy;
      result.theta = kin.theta * 180.0 / TMath::Pi();
      result.good = (result.kineticEnergy > 0.3 && result.kineticEnergy < 3.0 && result.theta > 70 && result.theta < 130);

      return result;
   }
};

TEST_F(ClusteringDigiScanTest, ScanDigitizedData)
{
   FairLogger::GetLogger()->SetLogScreenLevel("FATAL");

   // Try to open the digi file
   auto *file = TFile::Open("../macro/Simulation/ATTPC/16C_pp/data/output_digi.root");
   if (!file || file->IsZombie()) {
      std::cout << "SKIPPED: output_digi.root not found. Run C16_pp_sim.C + run_digi_attpc.C first." << std::endl;
      GTEST_SKIP();
   }

   auto *tree = (TTree *)file->Get("cbmsim");
   TClonesArray *patEvtArr = nullptr;
   tree->SetBranchAddress("AtPatternEvent", &patEvtArr);

   int nEvents = tree->GetEntries();

   struct Config {
      int method;       // 0=Smooth3D
      double radius;
      double distance;
      int straggling;   // 0=off, 1=on
      std::string label;
   };
   std::vector<Config> configs = {
      // Straggling OFF
      {0, 5.0, 15.0, 0, "r5d15 noStr"},
      {0, 10.0, 20.0, 0, "r10d20 noStr"},
      {0, 15.0, 20.0, 0, "r15d20 noStr"},
      {0, 15.0, 30.5, 0, "r15d30 noStr (def)"},
      {0, 20.0, 30.0, 0, "r20d30 noStr"},
      // Straggling ON (same configs)
      {0, 5.0, 15.0, 1, "r5d15 STR"},
      {0, 10.0, 20.0, 1, "r10d20 STR"},
      {0, 15.0, 20.0, 1, "r15d20 STR"},
      {0, 15.0, 30.5, 1, "r15d30 STR (def)"},
      {0, 20.0, 30.0, 1, "r20d30 STR"},
   };

   // Accumulate results per config
   struct ConfigStats {
      int nTried = 0;
      int nConverged = 0;
      int nGood = 0;
      double sumKE = 0;
      double sumClusters = 0;
   };
   std::vector<ConfigStats> stats(configs.size());

   // Event loop
   for (int iEvt = 0; iEvt < nEvents; iEvt++) {
      tree->GetEntry(iEvt);
      if (!patEvtArr || patEvtArr->GetEntries() == 0)
         continue;

      auto *patEvt = (AtPatternEvent *)patEvtArr->At(0);
      auto &tracks = patEvt->GetTrackCand();
      if (tracks.empty())
         continue;

      // Use largest track
      int bestTrack = 0;
      for (size_t t = 1; t < tracks.size(); t++) {
         if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
            bestTrack = t;
      }

      // Skip small tracks (beam events)
      if (tracks[bestTrack].GetHitArray().size() < 50)
         continue;

      // Scan clustering configs
      for (size_t c = 0; c < configs.size(); c++) {
         AtTrack trackCopy = tracks[bestTrack];
         trackCopy.ResetHitClusterArray();
         AtTools::AtTrackTransformer transformer;
         transformer.ClusterizeSmooth3D(trackCopy, configs[c].radius, configs[c].distance);
         auto result = FitTrack(trackCopy, iEvt, configs[c].straggling != 0);

         stats[c].nTried++;
         if (result.converged) {
            stats[c].nConverged++;
            stats[c].sumClusters += result.nClusters;
            if (result.good) {
               stats[c].nGood++;
               stats[c].sumKE += result.kineticEnergy;
            }
         }
      }
   }

   // Print summary table
   std::cout << "\n=============================================================" << std::endl;
   std::cout << " Clustering Scan on Digitized 16C(p,p) Data" << std::endl;
   std::cout << "=============================================================" << std::endl;
   std::cout << std::setw(20) << "config" << std::setw(8) << "tried" << std::setw(8) << "conv" << std::setw(8) << "good"
             << std::setw(10) << "good(%)" << std::setw(10) << "avgCl" << std::setw(10) << "avgKE" << std::endl;
   std::cout << std::string(74, '-') << std::endl;

   for (size_t c = 0; c < configs.size(); c++) {
      auto &s = stats[c];
      double goodPct = s.nTried > 0 ? 100.0 * s.nGood / s.nTried : 0;
      double avgCl = s.nConverged > 0 ? s.sumClusters / s.nConverged : 0;
      double avgKE = s.nGood > 0 ? s.sumKE / s.nGood : 0;
      std::cout << std::setw(20) << configs[c].label << std::setw(8) << s.nTried << std::setw(8) << s.nConverged
                << std::setw(8) << s.nGood << std::setw(10) << std::fixed << std::setprecision(1) << goodPct
                << std::setw(10) << std::setprecision(0) << avgCl << std::setw(10) << std::setprecision(2) << avgKE
                << std::endl;
   }
   std::cout << "=============================================================" << std::endl;

   // At least one config should give >50% good fits
   int bestGood = 0;
   for (auto &s : stats)
      bestGood = std::max(bestGood, s.nGood);
   EXPECT_GT(bestGood, 0) << "No clustering config produced any good fits";
}
