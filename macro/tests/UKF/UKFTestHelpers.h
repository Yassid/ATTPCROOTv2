/// @file UKFTestHelpers.h
/// @brief Shared helpers for UKF physics validation macros.
///
/// Include this file from each individual test macro:
///   #include "UKFTestHelpers.h"
///
/// Provides: constants, hit loading, UKF creation, and RunUKF().

#ifndef UKF_TEST_HELPERS_H
#define UKF_TEST_HELPERS_H

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using ROOT::Math::Polar3DVector;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const double mass_p = 938.272;           // Proton mass in MeV/c^2
const double charge_p = 1.602176634e-19; // Proton charge in C
const double gasDensity = 3.553e-5;      // g/cm^3 for 300 torr H2

// True initial state from GEANT simulation
const XYZPoint fTruePos(-3.40046e-04, -1.49863e-04, 1.0018);        // mm
const XYZVector fTrueMom(0.00935463e3, -0.0454279e3, 0.00826042e3); // MeV/c

// ---------------------------------------------------------------------------
// Hit loading
// ---------------------------------------------------------------------------
struct Hit {
   double x, y, z;
};

std::vector<Hit> gHits;

void LoadHits(int stride = 5)
{
   if (!gHits.empty())
      return;
   std::ifstream infile("hits.txt");
   if (!infile.is_open()) {
      std::cerr << "ERROR: Cannot open hits.txt" << std::endl;
      return;
   }
   double xi, yi, zi, ei;
   int count = 0;
   while (infile >> xi >> yi >> zi >> ei) {
      if (count % stride == 0)
         gHits.push_back({xi * 10, yi * 10, zi * 10}); // cm -> mm
      count++;
   }
   std::cout << "Loaded " << gHits.size() << " hits." << std::endl;
}

// ---------------------------------------------------------------------------
// UKF creation
// ---------------------------------------------------------------------------
kf::TrackFitterUKF *CreateUKF(double alpha = 1e-3, bool straggling = true)
{
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
   propagator.SetBField({0, 0, 2.85});

   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   auto *ukf = new kf::TrackFitterUKF(std::move(propagator), std::move(stepper));
   ukf->setParameters(alpha, 2.0, 0.0);
   ukf->fEnableEnStraggling = straggling;
   return ukf;
}

// ---------------------------------------------------------------------------
// Covariance helpers
// ---------------------------------------------------------------------------
TMatrixD MakeCov(double pMag)
{
   double sigma_pos = 1.0;
   double sigma_mom = 0.1 * pMag;
   double sigma_ang = 1.0 * M_PI / 180.0;
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 3; ++i)
      cov(i, i) = sigma_pos * sigma_pos;
   cov(3, 3) = sigma_mom * sigma_mom;
   cov(4, 4) = sigma_ang * sigma_ang;
   cov(5, 5) = sigma_ang * sigma_ang;
   return cov;
}

TMatrixD MakeMeasCov()
{
   TMatrixD cov(3, 3);
   cov.Zero();
   for (int i = 0; i < 3; ++i)
      cov(i, i) = 1.0;
   return cov;
}

// ---------------------------------------------------------------------------
// UKF runner
// ---------------------------------------------------------------------------
struct UKFResult {
   double pVertex;        // Smoothed vertex momentum
   double pVertexFilt;    // Filtered vertex momentum
   double chi2ndf;        // Chi2/ndf
   double rmsResidSmooth; // RMS smoothed residual
   double rmsResidFilt;   // RMS filtered residual
   int nTouch;            // Number of regularizations
   bool converged;
};

UKFResult RunUKF(kf::TrackFitterUKF *ukf, XYZPoint pos, XYZVector mom)
{
   UKFResult result{-1, -1, -1, -1, -1, 0, false};

   ukf->SetInitialState(pos, mom, MakeCov(mom.R()));
   ukf->SetMeasCov(MakeMeasCov());

   try {
      for (size_t i = 1; i < gHits.size(); ++i) {
         XYZPoint meas(gHits[i].x, gHits[i].y, gHits[i].z);
         ukf->predictUKF(meas);
         ukf->correctUKF(meas);
      }
      ukf->smoothUKF();
   } catch (const std::exception &e) {
      std::cout << "  Exception: " << e.what() << std::endl;
      return result;
   }

   result.converged = true;
   result.nTouch = ukf->nTouch;

   auto &smoothed = ukf->GetSmoothedStates();
   auto &filtered = ukf->GetFilteredStates();

   result.pVertex = smoothed[0][3];
   result.pVertexFilt = filtered[0][3];

   // Chi2/ndf
   double chi2 = 0;
   int n = smoothed.size();
   for (int i = 1; i < n; ++i) {
      XYZPoint sp(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
      XYZPoint mp(gHits[i].x, gHits[i].y, gHits[i].z);
      double d = (sp - mp).R();
      chi2 += d * d;
   }
   int ndf = std::max(1, 3 * (n - 1) - 6);
   result.chi2ndf = chi2 / ndf;

   // RMS residuals
   double sumFilt = 0, sumSmooth = 0;
   for (int i = 1; i < n; ++i) {
      XYZPoint mp(gHits[i].x, gHits[i].y, gHits[i].z);
      XYZPoint fp(filtered[i][0], filtered[i][1], filtered[i][2]);
      XYZPoint sp(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
      sumFilt += (fp - mp).Mag2();
      sumSmooth += (sp - mp).Mag2();
   }
   result.rmsResidFilt = std::sqrt(sumFilt / (n - 1));
   result.rmsResidSmooth = std::sqrt(sumSmooth / (n - 1));

   return result;
}

/// Call at the start of every test macro to suppress logging and load hits.
void InitTest()
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   LoadHits();
   if (gHits.size() < 10) {
      std::cerr << "Not enough hits loaded, aborting." << std::endl;
      gSystem->Exit(1);
   }
}

#endif // UKF_TEST_HELPERS_H
