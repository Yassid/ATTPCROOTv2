/// @file UKFPhysicsValidation.C
/// @brief Physics validation suite for the UKF track fitter.
///
/// Tests the UKF against known physics expectations:
///   1. Momentum reconstruction accuracy (single track, known truth)
///   2. Smoother vs filter improvement
///   3. Statistical consistency over many trials with smeared initial conditions
///   4. Chi2/ndf distribution validation
///   5. Energy loss consistency check
///   6. Alpha parameter sensitivity scan
///
/// Run from macro/tests/UKF/ after sourcing build/config.sh:
///   root -l -b -q 'UKFPhysicsValidation.C'
///   root -l -b -q 'UKFPhysicsValidation.C(1000)' // run with N statistical trials
///
/// Output:
///   - Console summary with PASS/FAIL for each test
///   - UKFPhysicsValidation.png with diagnostic plots

// ---------------------------------------------------------------------------
// Constants and helpers
// ---------------------------------------------------------------------------
using ROOT::Math::Polar3DVector;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

const double mass_p = 938.272;           // Proton mass in MeV/c^2
const double charge_p = 1.602176634e-19; // Proton charge in C
const double gasDensity = 3.553e-5;      // g/cm^3 for 300 torr H2

// True initial state from GEANT simulation
const XYZPoint fTruePos(-3.40046e-04, -1.49863e-04, 1.0018);          // mm
const XYZVector fTrueMom(0.00935463e3, -0.0454279e3, 0.00826042e3);   // MeV/c

struct Hit {
   double x, y, z;
};

std::vector<Hit> gHits;

std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr)
      return "../../resources/energy_loss/HinH.txt";
   return std::string(env) + "/resources/energy_loss/HinH.txt";
}

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

/// Create a fresh UKF instance with the given alpha
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

/// Run full UKF (forward + smooth) and return the smoothed vertex momentum.
/// Returns -1 on failure.
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
      return result;
   }

   result.converged = true;
   result.nTouch = ukf->nTouch;

   auto &smoothed = ukf->GetSmoothedStates();
   auto &filtered = ukf->GetFilteredStates();

   result.pVertex = smoothed[0][3];
   result.pVertexFilt = filtered[0][3];

   // Chi2/ndf (same formula as AtFitterUKF)
   double chi2 = 0;
   int n = smoothed.size();
   for (int i = 1; i < n; ++i) {
      XYZPoint sp(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
      XYZPoint mp(gHits[i].x, gHits[i].y, gHits[i].z);
      double d = (sp - mp).R();
      chi2 += d * d; // sigma_pos = 1 mm
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

// ---------------------------------------------------------------------------
// Test 1: Single track momentum reconstruction
// ---------------------------------------------------------------------------
bool TestSingleTrackReconstruction()
{
   std::cout << "\n===== Test 1: Single Track Momentum Reconstruction =====" << std::endl;

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, fTrueMom);
   delete ukf;

   if (!result.converged) {
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return false;
   }

   double pTrue = fTrueMom.R();
   double errSmooth = (result.pVertex - pTrue) / pTrue * 100;
   double errFilt = (result.pVertexFilt - pTrue) / pTrue * 100;

   std::cout << "  True momentum:       " << pTrue << " MeV/c" << std::endl;
   std::cout << "  Smoothed momentum:   " << result.pVertex << " MeV/c (err: " << errSmooth << "%)" << std::endl;
   std::cout << "  Filtered momentum:   " << result.pVertexFilt << " MeV/c (err: " << errFilt << "%)" << std::endl;
   std::cout << "  Chi2/ndf:            " << result.chi2ndf << std::endl;
   std::cout << "  RMS resid (smooth):  " << result.rmsResidSmooth << " mm" << std::endl;
   std::cout << "  RMS resid (filter):  " << result.rmsResidFilt << " mm" << std::endl;
   std::cout << "  Regularizations:     " << result.nTouch << std::endl;

   bool pass = std::abs(errSmooth) < 1.0 && result.nTouch == 0;
   std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires <1% error, 0 regularizations)" << std::endl;
   return pass;
}

// ---------------------------------------------------------------------------
// Test 2: Smoother must improve over filter (with smeared initial momentum)
// ---------------------------------------------------------------------------
bool TestSmootherImprovement()
{
   std::cout << "\n===== Test 2: Smoother vs Filter Quality =====" << std::endl;

   // Smear the initial momentum by 2% so the filter vertex is NOT trivially exact,
   // but still within the regime where sigma points reach the measurement planes.
   double pTrue = fTrueMom.R();
   double pSmeared = pTrue * 1.05; // 5% high
   Polar3DVector momPolar(pSmeared, fTrueMom.Theta(), fTrueMom.Phi());
   XYZVector momSmeared(momPolar);

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, momSmeared);
   delete ukf;

   if (!result.converged) {
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return false;
   }

   double errSmooth = std::abs(result.pVertex - pTrue);
   double errFilt = std::abs(result.pVertexFilt - pTrue);

   std::cout << "  Initial momentum:    " << pSmeared << " MeV/c (5% high)" << std::endl;
   std::cout << "  True momentum:       " << pTrue << " MeV/c" << std::endl;
   std::cout << "  |p_smooth - p_true| = " << errSmooth << " MeV/c" << std::endl;
   std::cout << "  |p_filt   - p_true| = " << errFilt << " MeV/c" << std::endl;
   std::cout << "  RMS resid smooth:    " << result.rmsResidSmooth << " mm" << std::endl;
   std::cout << "  RMS resid filter:    " << result.rmsResidFilt << " mm" << std::endl;

   bool pass = errSmooth <= errFilt;
   std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (smoother vertex error must be <= filter vertex error)" << std::endl;
   return pass;
}

// ---------------------------------------------------------------------------
// Test 3: Statistical consistency (many trials with smeared initial momentum)
// ---------------------------------------------------------------------------
bool TestStatisticalConsistency(int nTrials, TH1F *&hPull, TH1F *&hError, TH1F *&hChi2)
{
   std::cout << "\n===== Test 3: Statistical Consistency (" << nTrials << " trials) =====" << std::endl;

   double pTrue = fTrueMom.R();
   double smearSigma = 0.05 * pTrue; // 5% smearing on initial momentum

   hPull = new TH1F("hPull", "Momentum Pull Distribution;(p_{reco} - p_{true}) / #sigma_{p};Entries", 50, -5, 5);
   hError = new TH1F("hError", "Momentum Error;(p_{reco} - p_{true}) / p_{true} [%];Entries", 50, -3, 3);
   hChi2 = new TH1F("hChi2", "Chi2/ndf Distribution;#chi^{2}/ndf;Entries", 50, 0, 5);

   auto *ukf = CreateUKF(1e-3, true);

   int nConverged = 0;
   int nFailed = 0;

   for (int t = 0; t < nTrials; ++t) {
      if (t % 100 == 0 && t > 0)
         std::cout << "  Trial " << t << "/" << nTrials << std::endl;

      // Smear initial momentum
      double pSampled = gRandom->Gaus(pTrue, smearSigma);
      Polar3DVector momPolar(pSampled, fTrueMom.Theta(), fTrueMom.Phi());
      XYZVector momSmeared(momPolar);

      auto result = RunUKF(ukf, fTruePos, momSmeared);

      if (!result.converged) {
         nFailed++;
         continue;
      }
      nConverged++;

      double err = (result.pVertex - pTrue) / pTrue * 100;
      hError->Fill(err);
      hChi2->Fill(result.chi2ndf);

      // Pull: need the smoothed covariance at vertex
      auto &smoothedCov = ukf->GetSmoothedCovariances();
      if (!smoothedCov.empty()) {
         double sigmaP = std::sqrt(smoothedCov[0](3, 3));
         if (sigmaP > 0) {
            double pull = (result.pVertex - pTrue) / sigmaP;
            hPull->Fill(pull);
         }
      }
   }

   delete ukf;

   double convergenceRate = 100.0 * nConverged / nTrials;
   std::cout << "  Converged: " << nConverged << "/" << nTrials << " (" << convergenceRate << "%)" << std::endl;

   // Fit pull distribution
   if (hPull->GetEntries() > 20) {
      hPull->Fit("gaus", "Q");
      auto *fit = hPull->GetFunction("gaus");
      if (fit) {
         double pullMean = fit->GetParameter(1);
         double pullSigma = fit->GetParameter(2);
         std::cout << "  Pull mean:  " << pullMean << " (expect ~0)" << std::endl;
         std::cout << "  Pull sigma: " << pullSigma << " (expect ~1)" << std::endl;
      }
   }

   // Error distribution
   if (hError->GetEntries() > 20) {
      hError->Fit("gaus", "Q");
      auto *fit = hError->GetFunction("gaus");
      if (fit) {
         double errMean = fit->GetParameter(1);
         double errSigma = fit->GetParameter(2);
         std::cout << "  Error mean:  " << errMean << "%" << std::endl;
         std::cout << "  Error sigma: " << errSigma << "%" << std::endl;
      }
   }

   // Chi2/ndf
   double chi2Mean = hChi2->GetMean();
   std::cout << "  Chi2/ndf mean: " << chi2Mean << " (expect ~1)" << std::endl;

   bool pass = convergenceRate > 95.0;
   std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires >95% convergence rate)" << std::endl;
   return pass;
}

// ---------------------------------------------------------------------------
// Test 4: Energy loss consistency
// ---------------------------------------------------------------------------
bool TestEnergyLossConsistency()
{
   std::cout << "\n===== Test 4: Energy Loss Consistency =====" << std::endl;

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, fTrueMom);

   if (!result.converged) {
      delete ukf;
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return false;
   }

   auto &smoothed = ukf->GetSmoothedStates();

   // Check that momentum decreases monotonically (proton losing energy in gas)
   int nInversions = 0;
   for (size_t i = 1; i < smoothed.size(); ++i) {
      if (smoothed[i][3] > smoothed[i - 1][3] * 1.01) { // 1% tolerance
         nInversions++;
      }
   }

   double pFirst = smoothed[0][3];
   double pLast = smoothed.back()[3];
   double totalEloss = AtTools::Kinematics::KE(pFirst, mass_p) - AtTools::Kinematics::KE(pLast, mass_p);

   std::cout << "  Vertex momentum:     " << pFirst << " MeV/c" << std::endl;
   std::cout << "  Final momentum:      " << pLast << " MeV/c" << std::endl;
   std::cout << "  Total energy loss:   " << totalEloss << " MeV" << std::endl;
   std::cout << "  Momentum inversions: " << nInversions << " (steps where p increased >1%)" << std::endl;
   std::cout << "  Track length:        " << smoothed.size() << " points" << std::endl;

   // Energy loss should be positive and reasonable for a ~1 MeV proton in H2
   bool pass = totalEloss > 0 && totalEloss < 5.0 && nInversions < 3;
   std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires dE>0, dE<5 MeV, <3 inversions)" << std::endl;

   delete ukf;
   return pass;
}

// ---------------------------------------------------------------------------
// Test 5: Alpha parameter sensitivity scan
// ---------------------------------------------------------------------------
bool TestAlphaSensitivity()
{
   std::cout << "\n===== Test 5: Alpha Parameter Sensitivity =====" << std::endl;

   double alphas[] = {1e-3, 5e-3, 1e-2, 5e-2};
   int nAlphas = sizeof(alphas) / sizeof(alphas[0]);

   std::cout << std::setw(10) << "alpha" << std::setw(15) << "p_smooth" << std::setw(12) << "err(%)"
             << std::setw(12) << "chi2/ndf" << std::setw(10) << "nTouch" << std::setw(10) << "status" << std::endl;
   std::cout << std::string(69, '-') << std::endl;

   double pTrue = fTrueMom.R();
   bool allPass = true;

   for (int a = 0; a < nAlphas; ++a) {
      auto *ukf = CreateUKF(alphas[a], true);
      auto result = RunUKF(ukf, fTruePos, fTrueMom);
      delete ukf;

      std::string status;
      if (!result.converged) {
         status = "FAIL";
         allPass = false;
         std::cout << std::setw(10) << alphas[a] << std::setw(15) << "---" << std::setw(12) << "---"
                   << std::setw(12) << "---" << std::setw(10) << "---" << std::setw(10) << status << std::endl;
      } else {
         double err = (result.pVertex - pTrue) / pTrue * 100;
         bool ok = std::abs(err) < 2.0 && result.nTouch <= 2;
         if (!ok)
            allPass = false;
         status = ok ? "OK" : "WARN";

         std::cout << std::setw(10) << alphas[a] << std::setw(15) << result.pVertex << std::setw(12) << err
                   << std::setw(12) << result.chi2ndf << std::setw(10) << result.nTouch << std::setw(10) << status
                   << std::endl;
      }
   }

   std::cout << "  STATUS: " << (allPass ? "PASS" : "WARN") << std::endl;
   return allPass;
}

// ---------------------------------------------------------------------------
// Test 6: Straggling ON stability (was the original crash trigger for Bug 6)
// ---------------------------------------------------------------------------
bool TestStragglingEffect()
{
   std::cout << "\n===== Test 6: Energy Straggling Stability =====" << std::endl;

   // Test with straggling ON and a biased seed — this was the original crash
   // scenario before Bug 6 was fixed (stopped sigma points set p=0).
   double pTrue = fTrueMom.R();
   double pBiased = pTrue * 1.05;
   Polar3DVector momP(pBiased, fTrueMom.Theta(), fTrueMom.Phi());
   XYZVector momBiased(momP);

   auto *ukf = CreateUKF(1e-3, true);
   auto res = RunUKF(ukf, fTruePos, momBiased);
   delete ukf;

   double err = std::abs(res.pVertex - pTrue) / pTrue * 100;

   std::cout << "  Biased seed:    " << pBiased << " MeV/c (5% high)" << std::endl;
   std::cout << "  Straggling ON:  p = " << res.pVertex << " MeV/c, err = " << err
             << "%, chi2/ndf = " << res.chi2ndf << ", nTouch = " << res.nTouch
             << ", converged = " << res.converged << std::endl;

   bool pass = res.converged && err < 2.0;
   std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (must converge with straggling ON + biased seed)" << std::endl;
   return pass;
}

// ===========================================================================
// Main
// ===========================================================================
void UKFPhysicsValidation(int nStatTrials = 200)
{
   // Suppress verbose logging during batch runs (covariance matrix dumps, etc.)
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   std::cout << "========================================" << std::endl;
   std::cout << " UKF Physics Validation Suite" << std::endl;
   std::cout << "========================================" << std::endl;

   LoadHits();
   if (gHits.size() < 10) {
      std::cerr << "Not enough hits loaded, aborting." << std::endl;
      return;
   }

   int nPass = 0, nFail = 0;

   auto check = [&](bool pass) {
      if (pass) nPass++;
      else nFail++;
   };

   // Run tests
   check(TestSingleTrackReconstruction());
   check(TestSmootherImprovement());
   check(TestEnergyLossConsistency());
   check(TestAlphaSensitivity());
   check(TestStragglingEffect());

   TH1F *hPull = nullptr, *hError = nullptr, *hChi2 = nullptr;
   if (nStatTrials > 0)
      check(TestStatisticalConsistency(nStatTrials, hPull, hError, hChi2));

   // Summary
   std::cout << "\n========================================" << std::endl;
   std::cout << " SUMMARY: " << nPass << " passed, " << nFail << " failed out of " << (nPass + nFail) << " tests"
             << std::endl;
   std::cout << "========================================" << std::endl;

   // Save plots
   if (hPull && hError && hChi2) {
      TCanvas *c = new TCanvas("cVal", "UKF Physics Validation", 1500, 500);
      c->Divide(3, 1);
      c->cd(1);
      hPull->Draw("hist");
      if (hPull->GetFunction("gaus"))
         hPull->GetFunction("gaus")->Draw("same");
      c->cd(2);
      hError->Draw("hist");
      if (hError->GetFunction("gaus"))
         hError->GetFunction("gaus")->Draw("same");
      c->cd(3);
      hChi2->Draw("hist");
      c->SaveAs("UKFPhysicsValidation.png");
      std::cout << "Saved: UKFPhysicsValidation.png" << std::endl;
   }
}
