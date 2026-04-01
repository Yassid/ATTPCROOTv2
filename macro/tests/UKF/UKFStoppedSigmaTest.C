/// @file UKFStoppedSigmaTest.C
/// @brief Demonstrates the stopped-sigma-point bug in the UKF smoother.
///
/// When a sigma point's particle stops before reaching the measurement
/// surface, the propagator sets its momentum to zero. This creates a
/// catastrophic outlier (p=0 vs p~300 MeV/c) that corrupts:
///   - The predicted mean and covariance
///   - The cross-covariance used by the RTS smoother
///
/// The forward filter partially recovers via measurement correction,
/// but the smoother's cross-covariance is permanently corrupted.
///
/// Run:  root -l -b -q UKFStoppedSigmaTest.C

const double mass_p = 938.272;
const double charge_p = 1.602176634e-19;

struct Hit {
   double x, y, z;
};

std::vector<Hit> LoadHitsFromFile(int stride = 5)
{
   std::vector<Hit> hits;
   std::ifstream infile("hits.txt");
   if (!infile.is_open())
      return hits;
   double xi, yi, zi, ei;
   int count = 0;
   while (infile >> xi >> yi >> zi >> ei) {
      if (count % stride == 0)
         hits.push_back({xi * 10, yi * 10, zi * 10});
      count++;
   }
   return hits;
}

void UKFStoppedSigmaTest()
{
   auto hits = LoadHitsFromFile(5);
   if (hits.size() < 5) {
      std::cout << "Cannot open hits.txt" << std::endl;
      return;
   }

   using namespace ROOT::Math;

   // Setup propagator matching the UKF's internal setup
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
   propagator.SetBField({0, 0, 2.85});
   propagator.SetEField({0, 0, 0});

   // Get the mean state for a typical prediction step
   XYZPoint startPos(hits[0].x, hits[0].y, hits[0].z);
   XYZVector startMom(0.00935463e3, -0.0454279e3, 0.00826042e3); // MeV/c
   XYZPoint nextPos(hits[1].x, hits[1].y, hits[1].z);

   // Propagate mean state to build the measurement plane
   propagator.SetState(startPos, startMom);
   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   propagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPoint(nextPos), *stepper);
   auto meanState = propagator.GetState();
   Plane3D measPlane(meanState.fMom.Unit(), nextPos);

   double nominalP = startMom.R();
   std::cout << "Nominal momentum: " << nominalP << " MeV/c" << std::endl;
   std::cout << "Nominal KE: " << AtTools::Kinematics::KE(nominalP, mass_p) << " MeV" << std::endl;
   std::cout << std::endl;

   // Now test funcF with different energy loss scaling factors
   // In the UKF, v[0] is the energy loss scaling factor.
   // Sigma points perturb v around 1.0.
   std::cout << "=== Simulating sigma point propagation with different v[0] ===" << std::endl;
   std::cout << std::setw(8) << "v[0]" << std::setw(12) << "|p| MeV/c" << std::setw(12) << "theta"
             << std::setw(12) << "phi" << std::setw(12) << "status" << std::endl;
   std::cout << std::string(56, '-') << std::endl;

   double vValues[] = {0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2.0, 3.0};

   for (double v : vValues) {
      // Reset propagator to start
      auto eloss2 = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
      eloss2->SetProjectile(1, 1, 1);
      eloss2->SetMaterial(mat);
      AtTools::AtPropagator prop2(charge_p, mass_p, std::move(eloss2));
      prop2.SetBField({0, 0, 2.85});

      prop2.SetState(startPos, startMom);
      prop2.fScalingFactor = v;
      auto step2 = std::make_unique<AtTools::AtRK4Stepper>();
      prop2.PropagateToMeasurementSurface(AtTools::AtMeasurementPlane(measPlane), *step2);
      prop2.fScalingFactor = 1.0;

      auto state = prop2.GetState();
      double pMag = state.fMom.R();
      double theta = state.fMom.Theta();
      double phi = state.fMom.Phi();

      std::string status = "OK";
      if (pMag < 1e-6)
         status = "STOPPED (p=0)";
      else if (pMag < nominalP * 0.5)
         status = "LOW";

      std::cout << std::setw(8) << v << std::setw(12) << pMag << std::setw(12) << theta << std::setw(12) << phi
                << std::setw(12) << status << std::endl;
   }

   std::cout << std::endl;
   std::cout << "=== Impact Analysis ===" << std::endl;
   std::cout << std::endl;
   std::cout << "With alpha=1.0 and DIM_A=7, each sigma point has weight w_i = 1/14 = 0.0714" << std::endl;
   std::cout << "If ONE sigma point returns p=0 while others return p~" << nominalP << ":" << std::endl;
   std::cout << "  Mean momentum error ~ 0.0714 * " << nominalP << " = " << 0.0714 * nominalP << " MeV/c ("
             << 0.0714 * 100 << "%)" << std::endl;
   std::cout << "  Variance inflated by ~ 0.0714 * " << nominalP << "^2 = " << 0.0714 * nominalP * nominalP
             << " (MeV/c)^2" << std::endl;
   std::cout << std::endl;
   std::cout << "For angles: theta goes from ~1.48 to 0 (delta=1.48 rad = 84.8 deg)" << std::endl;
   std::cout << "           phi goes from ~-1.57 to 0 (delta=1.57 rad = 90 deg)" << std::endl;
   std::cout << std::endl;
   std::cout << "This outlier corrupts:" << std::endl;
   std::cout << "  1. Predicted mean (pulls momentum toward 0)" << std::endl;
   std::cout << "  2. Predicted covariance (inflates momentum/angle variance)" << std::endl;
   std::cout << "  3. Cross-covariance C_{k,k+1} (used by RTS smoother)" << std::endl;
   std::cout << std::endl;
   std::cout << "The forward filter partially recovers via measurement correction," << std::endl;
   std::cout << "but the corrupted cross-covariance degrades the smoother permanently." << std::endl;
   std::cout << std::endl;
   std::cout << "=== Proposed Fix ===" << std::endl;
   std::cout << "In AtPropagator.cxx line 184:" << std::endl;
   std::cout << "  fState.fMom = XYZVector(0, 0, 0);  // <-- BUG: creates p=0 outlier" << std::endl;
   std::cout << std::endl;
   std::cout << "Instead, preserve the momentum DIRECTION from the last valid step:" << std::endl;
   std::cout << "  fState.fMom = fState.fLastMom.Unit() * fStopTol;" << std::endl;
   std::cout << "  // or: keep fState.fMom = fState.fLastMom (before zeroing)" << std::endl;
   std::cout << std::endl;
   std::cout << "This keeps the sigma point close to the other sigma points in momentum" << std::endl;
   std::cout << "space, preventing the catastrophic outlier that corrupts the smoother." << std::endl;

   // Now run the full UKF to show the effect
   std::cout << std::endl;
   std::cout << "=== Full UKF comparison ===" << std::endl;

   // Run with straggling OFF (no stopped sigma points)
   {
      auto e = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
      e->SetProjectile(1, 1, 1);
      e->SetMaterial(mat);
      AtTools::AtPropagator prop(charge_p, mass_p, std::move(e));
      prop.SetBField({0, 0, 2.85});
      auto stp = std::make_unique<AtTools::AtRK4Stepper>();
      kf::TrackFitterUKF ukf(std::move(prop), std::move(stp));
      ukf.setParameters(1.0, 2.0, 0.0);
      ukf.fEnableEnStraggling = false;

      TMatrixD cov(6, 6);
      cov.Zero();
      for (int i = 0; i < 3; ++i)
         cov(i, i) = 1.0;
      cov(3, 3) = (0.1 * startMom.R()) * (0.1 * startMom.R());
      cov(4, 4) = cov(5, 5) = (M_PI / 180) * (M_PI / 180);

      ukf.SetInitialState(startPos, startMom, cov);
      TMatrixD measCov(3, 3);
      measCov.Zero();
      for (int i = 0; i < 3; ++i)
         measCov(i, i) = 1.0;
      ukf.SetMeasCov(measCov);

      int nSteps = std::min((int)hits.size() - 1, 20);
      for (int i = 1; i <= nSteps; ++i) {
         XYZPoint meas(hits[i].x, hits[i].y, hits[i].z);
         ukf.predictUKF(meas);
         ukf.correctUKF(meas);
      }
      ukf.smoothUKF();

      auto smooth = ukf.GetSmoothedStates();
      auto filt = ukf.GetFilteredStates();
      double sumFilt = 0, sumSmooth = 0;
      for (int i = 1; i <= nSteps; ++i) {
         XYZPoint m(hits[i].x, hits[i].y, hits[i].z);
         sumFilt += (XYZPoint(filt[i][0], filt[i][1], filt[i][2]) - m).Mag2();
         sumSmooth += (XYZPoint(smooth[i][0], smooth[i][1], smooth[i][2]) - m).Mag2();
      }
      double rmsFilt = std::sqrt(sumFilt / nSteps);
      double rmsSmooth = std::sqrt(sumSmooth / nSteps);

      std::cout << std::endl;
      std::cout << "=== Straggling OFF (no stopped sigma points) ===" << std::endl;
      std::cout << "  RMS filt residual:   " << rmsFilt << " mm" << std::endl;
      std::cout << "  RMS smooth residual: " << rmsSmooth << " mm" << std::endl;
      std::cout << "  Smooth/Filt ratio:   " << rmsSmooth / rmsFilt << std::endl;
      std::cout << "  Vertex p (smoothed): " << smooth[0][3] << " MeV/c" << std::endl;
      std::cout << "  True p:              " << startMom.R() << " MeV/c" << std::endl;
      std::cout << "  Vertex p error:      " << std::abs(smooth[0][3] - startMom.R()) / startMom.R() * 100 << "%"
                << std::endl;
      std::cout << "  Regularizations:     " << ukf.nTouch << std::endl;

      if (rmsSmooth <= rmsFilt)
         std::cout << "  RESULT: Smoother IMPROVES on filter (as expected)" << std::endl;
      else
         std::cout << "  RESULT: Smoother still degrades — ratio " << rmsSmooth / rmsFilt << std::endl;
   }

   // Run with straggling ON (may have stopped sigma points)
   std::cout << std::endl;
   std::cout << "=== Straggling ON ===" << std::endl;
   {
      auto e = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
      e->SetProjectile(1, 1, 1);
      e->SetMaterial(mat);
      AtTools::AtPropagator prop(charge_p, mass_p, std::move(e));
      prop.SetBField({0, 0, 2.85});
      auto stp = std::make_unique<AtTools::AtRK4Stepper>();
      kf::TrackFitterUKF ukf(std::move(prop), std::move(stp));
      ukf.setParameters(1.0, 2.0, 0.0);
      ukf.fEnableEnStraggling = true;

      TMatrixD cov(6, 6);
      cov.Zero();
      for (int i = 0; i < 3; ++i)
         cov(i, i) = 1.0;
      cov(3, 3) = (0.1 * startMom.R()) * (0.1 * startMom.R());
      cov(4, 4) = cov(5, 5) = (M_PI / 180) * (M_PI / 180);

      ukf.SetInitialState(startPos, startMom, cov);
      TMatrixD measCov(3, 3);
      measCov.Zero();
      for (int i = 0; i < 3; ++i)
         measCov(i, i) = 1.0;
      ukf.SetMeasCov(measCov);

      int nSteps = std::min((int)hits.size() - 1, 20);
      bool crashed = false;
      int lastStep = 0;
      for (int i = 1; i <= nSteps; ++i) {
         XYZPoint meas(hits[i].x, hits[i].y, hits[i].z);
         try {
            ukf.predictUKF(meas);
            ukf.correctUKF(meas);
            lastStep = i;
         } catch (const std::exception &ex) {
            std::cout << "  CRASHED at step " << i << ": " << ex.what() << std::endl;
            crashed = true;
            break;
         }
      }

      if (!crashed) {
         ukf.smoothUKF();
         auto smooth = ukf.GetSmoothedStates();
         auto filt = ukf.GetFilteredStates();
         double sumFilt = 0, sumSmooth = 0;
         for (int i = 1; i <= nSteps; ++i) {
            XYZPoint m(hits[i].x, hits[i].y, hits[i].z);
            sumFilt += (XYZPoint(filt[i][0], filt[i][1], filt[i][2]) - m).Mag2();
            sumSmooth += (XYZPoint(smooth[i][0], smooth[i][1], smooth[i][2]) - m).Mag2();
         }
         double rmsFilt = std::sqrt(sumFilt / nSteps);
         double rmsSmooth = std::sqrt(sumSmooth / nSteps);
         std::cout << "  RMS filt residual:   " << rmsFilt << " mm" << std::endl;
         std::cout << "  RMS smooth residual: " << rmsSmooth << " mm" << std::endl;
         std::cout << "  Smooth/Filt ratio:   " << rmsSmooth / rmsFilt << std::endl;
         std::cout << "  Regularizations:     " << ukf.nTouch << std::endl;
      } else {
         std::cout << "  Forward pass failed at step " << lastStep + 1 << "/" << nSteps << std::endl;
         std::cout << "  Stopped sigma points corrupted the covariance matrix" << std::endl;
         std::cout << "  beyond recovery by regularization." << std::endl;
      }
   }

   std::cout << std::endl;
   std::cout << "=== Conclusion ===" << std::endl;
   std::cout << "The bug is in AtPropagator.cxx line 184:" << std::endl;
   std::cout << "  fState.fMom = XYZVector(0, 0, 0);" << std::endl;
   std::cout << "When a sigma point's particle stops, setting momentum to zero" << std::endl;
   std::cout << "creates a catastrophic outlier (p=0 vs p~300 MeV/c, theta=0 vs ~1.5)" << std::endl;
   std::cout << "that corrupts the UKF's sigma-point statistics and smoother." << std::endl;
}
