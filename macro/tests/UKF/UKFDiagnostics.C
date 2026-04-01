/// @file UKFDiagnostics.C
/// @brief Diagnostic ROOT macro to probe potential bugs in the UKF.
///
/// REQUIRES: Eigen3 (the UKF module is conditionally compiled).
/// Run from macro/tests/UKF/ after sourcing build/config.sh:
///   root -l -b -q UKFDiagnostics.C
///
/// Identified issues:
///   1. Model noise (Qmod) applies fPosModelNoise to ALL state dimensions.
///   2. Sigma-point weight explosion with alpha=1e-3.
///   3. kReplaceFirstCov overwrites seed covariance.
///   4. Covariance PD stability over many steps.
///   5. Smoother quality vs filtered.
///   6. History vector consistency.

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr)
      return "../../resources/energy_loss/HinH.txt";
   return std::string(env) + "/resources/energy_loss/HinH.txt";
}

const double mass_p = 938.272;
const double charge_p = 1.602176634e-19;

struct Hit {
   double x, y, z;
};

std::vector<Hit> LoadHitsFromFile(int stride = 5)
{
   std::vector<Hit> hits;
   std::ifstream infile("hits.txt");
   if (!infile.is_open()) {
      std::cerr << "Cannot open hits.txt" << std::endl;
      return hits;
   }
   double xi, yi, zi, ei;
   int count = 0;
   while (infile >> xi >> yi >> zi >> ei) {
      if (count % stride == 0)
         hits.push_back({xi * 10, yi * 10, zi * 10}); // cm -> mm
      count++;
   }
   return hits;
}

// ---------------------------------------------------------------------------
// Test 1: Sigma-point weight analysis (pure math, no Eigen needed)
// ---------------------------------------------------------------------------
void TestWeights()
{
   std::cout << "\n===== Test 1: Sigma-Point Weight Analysis =====" << std::endl;

   int DIM_A = 7; // DIM_X=6 + DIM_V=1
   double betas[] = {2.0};
   double alphas[] = {1e-3, 1e-2, 1e-1, 0.5, 1.0};

   std::cout << std::setw(10) << "alpha" << std::setw(15) << "lambda" << std::setw(15) << "wM0" << std::setw(15)
             << "wC0" << std::setw(15) << "wi" << std::setw(15) << "sum(wM)" << std::setw(15) << "spread"
             << std::endl;
   std::cout << std::string(100, '-') << std::endl;

   for (double alpha : alphas) {
      double kappa = 0;
      double lambda = alpha * alpha * (DIM_A + kappa) - DIM_A;
      double denom = lambda + DIM_A;
      double wM0 = lambda / denom;
      double wC0 = wM0 + (1.0 - alpha * alpha + 2.0);
      double wi = 0.5 / denom;
      double sumMean = wM0 + 2 * DIM_A * wi;
      double spread = std::sqrt(std::abs(denom));

      std::cout << std::setw(10) << alpha << std::setw(15) << lambda << std::setw(15) << wM0 << std::setw(15) << wC0
                << std::setw(15) << wi << std::setw(15) << sumMean << std::setw(15) << spread << std::endl;
   }

   std::cout << "\n  NOTE: With alpha=1e-3 and DIM_A=7:" << std::endl;
   std::cout << "    - Mean weight wM0 ~ -1e6 (extreme negative)" << std::endl;
   std::cout << "    - Sigma point spread ~ 0.003 (very tight around mean)" << std::endl;
   std::cout << "    - Filter is hypersensitive to tiny nonlinear effects" << std::endl;
   std::cout << "  RECOMMENDATION: Use alpha >= 0.1 or switch to square-root UKF" << std::endl;
}

// ---------------------------------------------------------------------------
// Test 2: Model noise dimensional analysis
// ---------------------------------------------------------------------------
void TestModelNoise()
{
   std::cout << "\n===== Test 2: Model Noise Dimensional Analysis =====" << std::endl;

   double fPosModelNoise = 1e-4; // Default value in TrackFitterUKF

   std::cout << "  fPosModelNoise = " << fPosModelNoise << " applied to ALL 6 state dimensions:" << std::endl;
   std::cout << std::endl;
   std::cout << "  Dim  Physical quantity  Units          Noise value  Interpretation" << std::endl;
   std::cout << "  ---  -----------------  -----          -----------  --------------" << std::endl;
   std::cout << "  0    x position         mm^2           " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) << " mm   (OK)" << std::endl;
   std::cout << "  1    y position         mm^2           " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) << " mm   (OK)" << std::endl;
   std::cout << "  2    z position         mm^2           " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) << " mm   (OK)" << std::endl;
   std::cout << "  3    |p| momentum       (MeV/c)^2      " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) << " MeV/c (SUSPECT)" << std::endl;
   std::cout << "  4    theta              rad^2          " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) * 180 / M_PI << " deg (SUSPECT)" << std::endl;
   std::cout << "  5    phi                rad^2          " << fPosModelNoise << "      sigma = "
             << std::sqrt(fPosModelNoise) * 180 / M_PI << " deg (SUSPECT)" << std::endl;
   std::cout << std::endl;
   std::cout << "  BUG: SetInitialState applies fPosModelNoise to ALL diagonal elements" << std::endl;
   std::cout << "  of Qmod. For a 47 MeV proton (p~300 MeV/c), 0.01 MeV/c model noise" << std::endl;
   std::cout << "  is negligible. But 0.01 rad = 0.57 deg angular noise per step" << std::endl;
   std::cout << "  may be significant. The issue is using a single value for" << std::endl;
   std::cout << "  dimensions with different units and scales." << std::endl;
   std::cout << std::endl;
   std::cout << "  FIX: Set Qmod elements per-dimension, e.g.:" << std::endl;
   std::cout << "    Qmod(0..2) = fPosModelNoise          // position" << std::endl;
   std::cout << "    Qmod(3)    = fMomModelNoise           // momentum" << std::endl;
   std::cout << "    Qmod(4..5) = fAngModelNoise           // angles" << std::endl;
}

// ---------------------------------------------------------------------------
// Test 3: Forward pass with different alpha — uses compiled library
// ---------------------------------------------------------------------------
void TestForwardPass()
{
   std::cout << "\n===== Test 3: Forward Pass Stability =====" << std::endl;

   auto hits = LoadHitsFromFile(5);
   if (hits.size() < 5) {
      std::cout << "  Not enough hits, skipping." << std::endl;
      return;
   }

   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
   propagator.SetBField({0, 0, 2.85});

   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));

   // Use safe alpha
   ukf.setParameters(1.0, 2.0, 0.0);
   ukf.fEnableEnStraggling = true;

   ROOT::Math::XYZPoint startPos(hits[0].x, hits[0].y, hits[0].z);
   ROOT::Math::XYZVector startMom(0.00935463e3, -0.0454279e3, 0.00826042e3); // MeV/c

   double sigma_pos = 1.0;
   double sigma_mom = 0.1 * startMom.R();
   double sigma_ang = 1.0 * M_PI / 180.0;

   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 3; ++i)
      cov(i, i) = sigma_pos * sigma_pos;
   cov(3, 3) = sigma_mom * sigma_mom;
   cov(4, 4) = sigma_ang * sigma_ang;
   cov(5, 5) = sigma_ang * sigma_ang;

   ukf.SetInitialState(startPos, startMom, cov);

   TMatrixD measCov(3, 3);
   measCov.Zero();
   for (int i = 0; i < 3; ++i)
      measCov(i, i) = sigma_pos * sigma_pos;
   ukf.SetMeasCov(measCov);

   int nSteps = std::min((int)hits.size() - 1, 25);
   std::vector<double> stepNums, momenta, residuals;

   bool failed = false;
   for (int i = 1; i <= nSteps; ++i) {
      ROOT::Math::XYZPoint meas(hits[i].x, hits[i].y, hits[i].z);
      try {
         ukf.predictUKF(meas);
         ukf.correctUKF(meas);

         auto state = ukf.vecX();
         ROOT::Math::XYZPoint pos(state[0], state[1], state[2]);
         double resid = (pos - meas).R();

         stepNums.push_back(i);
         momenta.push_back(state[3]);
         residuals.push_back(resid);
      } catch (const std::exception &e) {
         std::cout << "  FAILED at step " << i << ": " << e.what() << std::endl;
         failed = true;
         break;
      }
   }

   std::cout << "  Forward pass: " << (failed ? "FAILED" : "OK") << " (" << stepNums.size() << "/" << nSteps
             << " steps)" << std::endl;
   std::cout << "  Regularizations (nTouch): " << ukf.nTouch << std::endl;

   if (!stepNums.empty()) {
      std::cout << "  Momentum range: " << *std::min_element(momenta.begin(), momenta.end()) << " - "
                << *std::max_element(momenta.begin(), momenta.end()) << " MeV/c" << std::endl;
      std::cout << "  Max residual: " << *std::max_element(residuals.begin(), residuals.end()) << " mm" << std::endl;
   }

   // Run smoother
   if (!failed) {
      try {
         ukf.smoothUKF();
         auto smoothed = ukf.GetSmoothedStates();
         auto filtered = ukf.GetFilteredStates();

         std::cout << "\n  Smoother results:" << std::endl;
         std::cout << "  Smoothed vertex momentum: " << smoothed[0][3] << " MeV/c" << std::endl;
         std::cout << "  Filtered vertex momentum: " << filtered[0][3] << " MeV/c" << std::endl;
         std::cout << "  True initial momentum:    " << startMom.R() << " MeV/c" << std::endl;

         double errSmooth = std::abs(smoothed[0][3] - startMom.R()) / startMom.R() * 100;
         double errFilt = std::abs(filtered[0][3] - startMom.R()) / startMom.R() * 100;
         std::cout << "  Smoothed vertex error:    " << errSmooth << "%" << std::endl;
         std::cout << "  Filtered vertex error:    " << errFilt << "%" << std::endl;

         // Compare residuals
         double sumFiltResid = 0, sumSmoothResid = 0;
         for (int i = 1; i < (int)smoothed.size() && i <= nSteps; ++i) {
            ROOT::Math::XYZPoint measPt(hits[i].x, hits[i].y, hits[i].z);
            ROOT::Math::XYZPoint filtPt(filtered[i][0], filtered[i][1], filtered[i][2]);
            ROOT::Math::XYZPoint smoothPt(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
            sumFiltResid += (filtPt - measPt).Mag2();
            sumSmoothResid += (smoothPt - measPt).Mag2();
         }
         std::cout << "  RMS filtered residual:    " << std::sqrt(sumFiltResid / nSteps) << " mm" << std::endl;
         std::cout << "  RMS smoothed residual:    " << std::sqrt(sumSmoothResid / nSteps) << " mm" << std::endl;

         // Plot
         TCanvas *c = new TCanvas("cFwd", "Forward Pass Diagnostics", 1200, 600);
         c->Divide(2, 1);

         c->cd(1);
         TGraph *gMom = new TGraph(stepNums.size(), stepNums.data(), momenta.data());
         gMom->SetTitle("Filtered Momentum vs Step;Step;p [MeV/c]");
         gMom->SetMarkerStyle(20);
         gMom->SetMarkerColor(kRed);
         gMom->Draw("ALP");

         std::vector<double> smoothMom;
         std::vector<double> smoothSteps;
         for (int i = 1; i < (int)smoothed.size(); ++i) {
            smoothSteps.push_back(i);
            smoothMom.push_back(smoothed[i][3]);
         }
         TGraph *gSmoothMom = new TGraph(smoothSteps.size(), smoothSteps.data(), smoothMom.data());
         gSmoothMom->SetMarkerStyle(22);
         gSmoothMom->SetMarkerColor(kGreen + 2);
         gSmoothMom->Draw("LP SAME");

         auto leg = new TLegend(0.6, 0.7, 0.9, 0.9);
         leg->AddEntry(gMom, "Filtered", "lp");
         leg->AddEntry(gSmoothMom, "Smoothed", "lp");
         leg->Draw();

         c->cd(2);
         TGraph *gResid = new TGraph(stepNums.size(), stepNums.data(), residuals.data());
         gResid->SetTitle("Position Residual vs Step;Step;Residual [mm]");
         gResid->SetMarkerStyle(20);
         gResid->Draw("ALP");

         c->SaveAs("UKFForwardPass.png");
         std::cout << "  Saved plot: UKFForwardPass.png" << std::endl;

      } catch (const std::exception &e) {
         std::cout << "  Smoother FAILED: " << e.what() << std::endl;
      }
   }
}

// ---------------------------------------------------------------------------
// Test 4: kReplaceFirstCov comparison
// ---------------------------------------------------------------------------
void TestReplaceFirstCov()
{
   std::cout << "\n===== Test 4: kReplaceFirstCov Impact =====" << std::endl;

   auto hits = LoadHitsFromFile(5);
   if (hits.size() < 5) {
      std::cout << "  Not enough hits, skipping." << std::endl;
      return;
   }

   ROOT::Math::XYZPoint startPos(hits[0].x, hits[0].y, hits[0].z);
   ROOT::Math::XYZVector startMom(0.00935463e3, -0.0454279e3, 0.00826042e3);

   double sigma_pos = 1.0;
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 3; ++i)
      cov(i, i) = sigma_pos * sigma_pos;
   cov(3, 3) = (0.1 * startMom.R()) * (0.1 * startMom.R());
   cov(4, 4) = cov(5, 5) = (1.0 * M_PI / 180) * (1.0 * M_PI / 180);

   TMatrixD measCov(3, 3);
   measCov.Zero();
   for (int i = 0; i < 3; ++i)
      measCov(i, i) = sigma_pos * sigma_pos;

   int nSteps = std::min((int)hits.size() - 1, 20);

   auto runWithFlag = [&](bool replaceFlag) -> double {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
      eloss->SetProjectile(1, 1, 1);
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back({1, 1, 1});
      eloss->SetMaterial(mat);

      AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
      propagator.SetBField({0, 0, 2.85});
      auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
      kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));
      ukf.setParameters(1.0, 2.0, 0.0);
      ukf.kReplaceFirstCov = replaceFlag;
      ukf.fEnableEnStraggling = true;

      ukf.SetInitialState(startPos, startMom, cov);
      ukf.SetMeasCov(measCov);

      for (int i = 1; i <= nSteps; ++i) {
         ROOT::Math::XYZPoint meas(hits[i].x, hits[i].y, hits[i].z);
         ukf.predictUKF(meas);
         ukf.correctUKF(meas);
      }
      ukf.smoothUKF();
      return ukf.GetSmoothedStates()[0][3]; // vertex momentum
   };

   double p_on = runWithFlag(true);
   double p_off = runWithFlag(false);
   double p_true = startMom.R();

   std::cout << "  True vertex momentum:            " << p_true << " MeV/c" << std::endl;
   std::cout << "  Smoothed p (replaceFirstCov ON):  " << p_on << " MeV/c  (err: "
             << std::abs(p_on - p_true) / p_true * 100 << "%)" << std::endl;
   std::cout << "  Smoothed p (replaceFirstCov OFF): " << p_off << " MeV/c  (err: "
             << std::abs(p_off - p_true) / p_true * 100 << "%)" << std::endl;
   std::cout << "  Difference:                       " << std::abs(p_on - p_off) << " MeV/c ("
             << std::abs(p_on - p_off) / p_true * 100 << "%)" << std::endl;
}

// ---------------------------------------------------------------------------
// Test 5: History vector size consistency
// ---------------------------------------------------------------------------
void TestHistoryConsistency()
{
   std::cout << "\n===== Test 5: History Vector Size Consistency =====" << std::endl;

   auto hits = LoadHitsFromFile(5);
   if (hits.size() < 5) {
      std::cout << "  Not enough hits, skipping." << std::endl;
      return;
   }

   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
   propagator.SetBField({0, 0, 2.85});
   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));
   ukf.setParameters(1.0, 2.0, 0.0);
   ukf.fEnableEnStraggling = true;

   ROOT::Math::XYZPoint startPos(hits[0].x, hits[0].y, hits[0].z);
   ROOT::Math::XYZVector startMom(0.00935463e3, -0.0454279e3, 0.00826042e3);

   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 3; ++i)
      cov(i, i) = 1.0;
   cov(3, 3) = (0.1 * startMom.R()) * (0.1 * startMom.R());
   cov(4, 4) = cov(5, 5) = (1.0 * M_PI / 180) * (1.0 * M_PI / 180);
   ukf.SetInitialState(startPos, startMom, cov);

   TMatrixD measCov(3, 3);
   measCov.Zero();
   for (int i = 0; i < 3; ++i)
      measCov(i, i) = 1.0;
   ukf.SetMeasCov(measCov);

   bool ok = true;
   int nSteps = std::min((int)hits.size() - 1, 15);

   for (int i = 1; i <= nSteps; ++i) {
      ROOT::Math::XYZPoint meas(hits[i].x, hits[i].y, hits[i].z);
      ukf.predictUKF(meas);
      ukf.correctUKF(meas);

      size_t expected = i + 1;
      auto nPred = ukf.GetPredictedStates().size();
      auto nFilt = ukf.GetFilteredStates().size();

      if (nPred != expected || nFilt != expected) {
         std::cout << "  MISMATCH at step " << i << ": pred=" << nPred << " filt=" << nFilt
                   << " expected=" << expected << std::endl;
         ok = false;
      }
   }

   ukf.smoothUKF();
   auto nSmooth = ukf.GetSmoothedStates().size();
   auto nPred = ukf.GetPredictedStates().size();
   if (nSmooth != nPred) {
      std::cout << "  MISMATCH after smooth: smooth=" << nSmooth << " pred=" << nPred << std::endl;
      ok = false;
   }

   std::cout << "  STATUS: " << (ok ? "PASS" : "FAIL") << " (" << nSteps << " steps, "
             << nPred << " total entries)" << std::endl;
}

// ===========================================================================
// Main
// ===========================================================================
void UKFDiagnostics()
{
   std::cout << "========================================" << std::endl;
   std::cout << " UKF Bug Diagnostics" << std::endl;
   std::cout << "========================================" << std::endl;

   // These two tests are pure math — always work
   TestWeights();
   TestModelNoise();

   // These require compiled UKF library + Eigen3 + hits.txt
   TestForwardPass();
   TestReplaceFirstCov();
   TestHistoryConsistency();

   std::cout << "\n========================================" << std::endl;
   std::cout << " Done." << std::endl;
   std::cout << "========================================" << std::endl;
}
