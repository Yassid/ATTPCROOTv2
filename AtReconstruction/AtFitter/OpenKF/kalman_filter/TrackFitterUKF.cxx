#include "TrackFitterUKF.h"

#include <FairLogger.h>

#include <Math/Plane3D.h>
#include <TMatrixD.h>

#include <AtPropagator.h>

namespace kf {

void TrackFitterUKF::Reset()
{
   // Reset the state vector and covariance matrix
   TrackFitterUKFBase::Reset();
   m_vecXSmooth.clear();
   m_matPSmooth.clear();
   m_matCPredHist.clear();
   fMeanStep = AtTools::AtPropagator::StepState(); // Reset the step state
}

void TrackFitterUKF::SetInitialState(const ROOT::Math::XYZPoint &initialPosition,
                                     const ROOT::Math::XYZVector &initialMomentum, const TMatrixD &initialCovariance)
{
   // If we are setting the initial state, then we should clear the history.
   Reset();
   fPropagator.SetState(initialPosition, initialMomentum); // Set the initial state in the propagator
   m_vecX[0] = initialPosition.X();                        // X position
   m_vecX[1] = initialPosition.Y();                        // Y position
   m_vecX[2] = initialPosition.Z();                        // Z position
   m_vecX[3] = initialMomentum.R();                        // Momentum magnitude
   m_vecX[4] = initialMomentum.Theta();                    // Polar angle
   m_vecX[5] = initialMomentum.Phi();                      // Azimuthal angle

   // Copy elements from initialCovariance to m_matP
   for (int i = 0; i < m_matP.rows(); ++i) {
      for (int j = 0; j < m_matP.cols(); ++j) {
         m_matP(i, j) = initialCovariance(i, j);
      }
   }

   // Set model noise per-dimension: position (mm^2), momentum (MeV/c)^2, angles (rad^2)
   m_matQmod(0, 0) = fPosModelNoise; // X position
   m_matQmod(1, 1) = fPosModelNoise; // Y position
   m_matQmod(2, 2) = fPosModelNoise; // Z position
   m_matQmod(3, 3) = fMomModelNoise; // Momentum magnitude
   m_matQmod(4, 4) = fAngModelNoise; // Theta
   m_matQmod(5, 5) = fAngModelNoise; // Phi

   // Save the initial state in our history vectors
   m_vecXFiltHist.push_back(m_vecX);
   m_matPFiltHist.push_back(m_matP);
   m_vecXPredHist.push_back(m_vecX);
   m_matPPredHist.push_back(m_matP);
   m_matCPredHist.push_back(Matrix<TF_DIM_X, TF_DIM_X>::Zero()); // Cross-correlation is not defined for the first point

   // We need to calculate the sigma points for the initial state
   updateAugmentedStateAndCovariance();                   // Update the augmented state vector and covariance matrix
   m_matSigmaXa = calculateSigmaPoints(m_vecXa, m_matPa); // Calculate the sigma points for the initial state
   // Now we grab the sigma points only for the state.
   m_matSigmaXPred = m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A); // Extract the state sigma points

   logEigen("Initial cov", m_matP, 0); // Log the eigenvalues of the initial covariance matrix
   LOG(info) << "Initial COV:" << std::endl << m_matP;
}

TMatrixD TrackFitterUKF::GetStateCovariance() const
{
   TMatrixD cov(m_matP.rows(), m_matP.cols());
   for (int i = 0; i < m_matP.rows(); ++i) {
      for (int j = 0; j < m_matP.cols(); ++j) {
         cov(i, j) = m_matP(i, j); // Copy covariance matrix to TMatrixD
      }
   }
   return cov;
}

std::array<double, TrackFitterUKF::DIM_A> TrackFitterUKF::GetAugStateVector() const
{
   std::array<double, DIM_A> result;
   for (int i = 0; i < DIM_A; ++i) {
      result[i] = m_vecXa[i];
   }
   return result;
}
TMatrixD TrackFitterUKF::GetAugStateCovariance() const
{
   TMatrixD cov(m_matPa.rows(), m_matPa.cols());
   for (int i = 0; i < m_matPa.rows(); ++i) {
      for (int j = 0; j < m_matPa.cols(); ++j) {
         cov(i, j) = m_matPa(i, j); // Copy augmented covariance matrix to TMatrixD
      }
   }
   return cov;
}

void TrackFitterUKF::SetMeasCov(const TMatrixD &measCov)
{
   if (measCov.GetNrows() != TF_DIM_Z || measCov.GetNcols() != TF_DIM_Z) {
      throw std::runtime_error("Measurement covariance matrix must be of size " + std::to_string(TF_DIM_Z) + "x" +
                               std::to_string(TF_DIM_Z));
   }
   for (int i = 0; i < TF_DIM_Z; ++i) {
      for (int j = 0; j < TF_DIM_Z; ++j) {
         m_matR(i, j) = measCov(i, j); // Copy measurement covariance to m_matR
      }
   }
}

std::array<float32_t, TrackFitterUKF::TF_DIM_V> TrackFitterUKF::calculateProcessNoiseMean()
{
   // The process noise is the scaling factor for dedx. By definition the mean should be 1
   std::array<float32_t, TF_DIM_V> processNoiseMean{1};
   return processNoiseMean;
}

Matrix<TrackFitterUKF::TF_DIM_V, TrackFitterUKF::TF_DIM_V> TrackFitterUKF::calculateProcessNoiseCovariance()
{
   assert(TF_DIM_V == 1 && "Process noise covariance is only implemented for DIM_V = 1");
   // Calculate the process noise covariance matrix
   Matrix<TF_DIM_V, TF_DIM_V> matQ{Matrix<TF_DIM_V, TF_DIM_V>::Zero()};

   // We need to know what the energy of the particle before/after transport.
   double eIn = AtTools::Kinematics::KE(fMeanStep.fLastMom, fMeanStep.fMass);
   double eOut = AtTools::Kinematics::KE(fMeanStep.fMom, fMeanStep.fMass);

   if (!fEnableEnStraggling) {
      // Set a small floor rather than exactly zero. A zero variance makes the
      // augmented covariance matrix singular, which causes the Cholesky
      // decomposition in calculateSigmaPoints() to fail.
      matQ(0, 0) = 1e-12;
      return matQ;
   }

   if (const auto *elossModel = fPropagator.GetELossModel()) {
      double dedx_straggle = elossModel->GetdEdxStraggling(eIn, eOut);
      double factor = dedx_straggle / elossModel->GetdEdx(eIn);
      if (factor > fMaxStragglingFactor) {
         LOG(warn) << "Process noise factor for energy straggling is greater than " << fMaxStragglingFactor
                   << ". To maintain stability, we will "
                      "use a factor of "
                   << fMaxStragglingFactor << ".";
         factor = fMaxStragglingFactor;
      }
      matQ(0, 0) = factor * factor; // Variance for the dedx straggling.
      LOG(debug) << "Calculating process noise for straggling between " << eIn << " MeV and " << eOut << " MeV over "
                 << elossModel->GetRange(eIn, eOut) << " mm.";
      LOG(debug) << "Process noise covariance for energy straggling: " << matQ(0, 0) << " (factor: " << factor
                 << ", dedx_straggle: " << dedx_straggle << ", dEdx: " << elossModel->GetdEdx(eIn) << ")";

   } else {
      throw std::runtime_error("Cannot calculate process noise covariance without an energy loss model");
   }
   // TODO: Add multiple scattering
   return matQ;
}

Vector<TrackFitterUKF::TF_DIM_X> TrackFitterUKF::funcF(const Vector<TrackFitterUKF::TF_DIM_X> &x,
                                                       const Vector<TrackFitterUKF::TF_DIM_V> &v,
                                                       const Vector<TrackFitterUKF::TF_DIM_Z> &z)
{
   // Set the state of the propagator to the current state vector
   using namespace ROOT::Math;

   // Validate sigma point state: momentum must be positive and angles finite
   double pMag = x[3];
   double theta = x[4];
   double phi = x[5];
   if (pMag <= 0 || std::isnan(pMag) || std::isinf(pMag) || std::isnan(theta) || std::isnan(phi)) {
      // Return the input unchanged — this sigma point is invalid
      Vector<TF_DIM_X> vecX;
      for (int d = 0; d < TF_DIM_X; d++)
         vecX[d] = x[d];
      return vecX;
   }

   XYZPoint fPos(x[0], x[1], x[2]);
   Polar3DVector fMom(pMag, theta, phi);

   fPropagator.SetState(fPos, XYZVector(fMom));

   fPropagator.fScalingFactor = v[0] * fELossScaleFactor;
   fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPlane(fMeasurementPlane), *fStepper);
   fPropagator.fScalingFactor = fELossScaleFactor;

   auto fState = fPropagator.GetState();
   Vector<TF_DIM_X> vecX{Vector<TF_DIM_X>::Zero()};
   vecX[0] = fState.fPos.X();
   vecX[1] = fState.fPos.Y();
   vecX[2] = fState.fPos.Z();
   vecX[3] = fState.fMom.R();
   vecX[4] = fState.fMom.Theta();
   vecX[5] = fState.fMom.Phi();

   // Validate output
   for (int d = 0; d < TF_DIM_X; d++) {
      if (std::isnan(vecX[d]) || std::isinf(vecX[d]))
         vecX[d] = x[d]; // Fall back to input
   }

   return vecX;
}

Vector<TrackFitterUKF::TF_DIM_Z> TrackFitterUKF::funcH(const Vector<TrackFitterUKF::TF_DIM_X> &x)
{
   Vector<TF_DIM_Z> vecZ;
   // Calculate the measurement vector based on the position and momentum
   vecZ[0] = x[0]; // X coordinate
   vecZ[1] = x[1]; // Y coordinate
   vecZ[2] = x[2]; // Z coordinate

   return vecZ; // Return the measurement vector
}

void TrackFitterUKF::predictUKF(const ROOT::Math::XYZPoint &z)
{
   using namespace ROOT::Math;

   // First we need to propagate the mean state vector to the next measurement point.
   XYZPoint startingPosition{m_vecX[0], m_vecX[1], m_vecX[2]};      // Get the starting position from the state vector
   Polar3DVector startingMomentum{m_vecX[3], m_vecX[4], m_vecX[5]}; // Get the starting momentum from the state vector

   LOG(debug) << "Propagating reference state from position: " << startingPosition
              << " with momentum: " << XYZVector(startingMomentum);

   fPropagator.SetState(startingPosition, XYZVector(startingMomentum));
   fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPoint(z), *fStepper);
   fMeanStep = fPropagator.GetState();    // Get the mean step information from the propagator
   fMeanStep.fLastPos = startingPosition; // Store the last position
   fMeanStep.fLastMom = startingMomentum; // Store the last momentum

   LOG(debug) << "Propagated to position: " << fMeanStep.fPos << " with momentum: " << fMeanStep.fMom;

   // Now we can construct the reference plane.
   fMeasurementPlane = Plane3D(fMeanStep.fMom.Unit(),
                               XYZPoint(z)); // Create a plane using the momentum direction and position
   Vector<TF_DIM_Z> zVec;                    // Initialize the measurement vector
   zVec[0] = z.X();
   zVec[1] = z.Y();
   zVec[2] = z.Z();
   auto callback = [this](const kf::Vector<TF_DIM_X> &x_, const kf::Vector<TF_DIM_V> &v_,
                          const kf::Vector<TF_DIM_Z> &z_) { return funcF(x_, v_, z_); };

   TrackFitterUKFBase::predictUKF(callback, zVec);

   // Get the sigma points belonging to the predicted state
   Matrix<TF_DIM_X, SIGMA_DIM_A> sigmaXx{m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A)};

   // Calculate the cross-corelation between the filtered state at k and predicted state at k+1
   auto matCPred =
      calculateCrossCorrelation<TF_DIM_X>(m_matSigmaXPred, m_vecXFiltHist.back(), sigmaXx, m_vecXPredHist.back());
   m_matCPredHist.push_back(matCPred); // Store the cross-correlation matrix
}

void TrackFitterUKF::predictUKFRef(const ROOT::Math::XYZPoint &z, const std::array<double, TF_DIM_X> &xRefPrev)
{
   using namespace ROOT::Math;

   // Current filtered estimate at i-1 (its deviation from the reference is what we transport).
   Vector<TF_DIM_X> xEst = m_vecX;
   Vector<TF_DIM_X> xRefVec;
   for (int d = 0; d < TF_DIM_X; ++d)
      xRefVec[d] = xRefPrev[d];

   // Measurement plane from propagating the REFERENCE state to z (linearization surface).
   XYZPoint refPos{xRefPrev[0], xRefPrev[1], xRefPrev[2]};
   Polar3DVector refMom{xRefPrev[3], xRefPrev[4], xRefPrev[5]};
   fPropagator.SetState(refPos, XYZVector(refMom));
   fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPoint(z), *fStepper);
   auto refStep = fPropagator.GetState();
   fMeasurementPlane = Plane3D(refStep.fMom.Unit(), XYZPoint(z));

   Vector<TF_DIM_Z> zVec;
   zVec[0] = z.X();
   zVec[1] = z.Y();
   zVec[2] = z.Z();

   // Augmented state/cov: state block CENTERED ON THE REFERENCE (spread = current cov), noise block.
   m_vecXa.setZero();
   for (int d = 0; d < TF_DIM_X; ++d)
      m_vecXa[d] = xRefVec[d];
   auto pnMean = calculateProcessNoiseMean();
   for (int d = 0; d < TF_DIM_V; ++d)
      m_vecXa[TF_DIM_X + d] = pnMean[d];
   m_matPa.setZero();
   for (int a = 0; a < TF_DIM_X; ++a)
      for (int b = 0; b < TF_DIM_X; ++b)
         m_matPa(a, b) = m_matP(a, b);
   m_matQaug = calculateProcessNoiseCovariance();
   for (int a = 0; a < TF_DIM_V; ++a)
      for (int b = 0; b < TF_DIM_V; ++b)
         m_matPa(TF_DIM_X + a, TF_DIM_X + b) = m_matQaug(a, b);

   m_matSigmaXa = calculateSigmaPoints(m_vecXa, m_matPa);
   Matrix<TF_DIM_X, SIGMA_DIM_A> sigmaXin{m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A)};
   Matrix<TF_DIM_V, SIGMA_DIM_A> sigmaXv{m_matSigmaXa.block(TF_DIM_X, 0, TF_DIM_V, SIGMA_DIM_A)};

   // Propagate each state sigma point (linearized around the reference).
   Matrix<TF_DIM_X, SIGMA_DIM_A> sigmaY;
   for (int32_t i = 0; i < SIGMA_DIM_A; ++i) {
      const Vector<TF_DIM_X> xi{sigmaXin.col(i)};
      const Vector<TF_DIM_V> vi{sigmaXv.col(i)};
      sigmaY.col(i) = funcF(xi, vi, zVec);
   }

   // Predicted reference mean + covariance from propagated sigma points.
   Vector<TF_DIM_X> yRef;
   Matrix<TF_DIM_X, TF_DIM_X> Pyy;
   calculateWeightedMeanAndCovariance<TF_DIM_X>(sigmaY, yRef, Pyy);

   // Cross-corr (for the smoother) from the reference-centered sigma: Pxy = Cov(Xin, Y).
   Matrix<TF_DIM_X, TF_DIM_X> Pxy = calculateCrossCorrelation<TF_DIM_X>(sigmaXin, xRefVec, sigmaY, yRef);

   // Predicted MEAN = full nonlinear propagation of the actual estimate (funcF(xEst)) —
   // as accurate as the standard predict, with NO linearization error. The earlier
   // yRef + F*(xEst-xRef) linearization underestimated the nonlinear propagation and,
   // re-applied each iteration, compounded into divergence (0.654->0.725->0.826).
   // The reference is used only to STABILIZE the covariance/cross-corr (sigma stay in
   // the linear regime), not the mean. shift recenters the sigma onto the mean.
   Vector<TF_DIM_V> pnVec;
   for (int d = 0; d < TF_DIM_V; ++d)
      pnVec[d] = pnMean[d];
   Vector<TF_DIM_X> yEst = funcF(xEst, pnVec, zVec);
   bool badMean = false;
   for (int d = 0; d < TF_DIM_X; ++d)
      if (!std::isfinite(yEst[d])) badMean = true;
   if (yEst[3] <= 0 || yEst[3] > 1000.0) badMean = true;
   Vector<TF_DIM_X> mean = badMean ? yRef : yEst; // fall back to reference propagation if estimate prop is bad
   Vector<TF_DIM_X> shift = mean - yRef;

   if (std::getenv("REFDBG")) {
      double dn = (xEst - xRefVec).norm();
      fprintf(stderr, "[REF] |d|=%.3g refp=%.3g yRefp=%.3g meanp=%.3g badMean=%d Pyy33=%.3g\n", dn, xRefVec[3], yRef[3],
              mean[3], (int)badMean, Pyy(3, 3));
   }

   // Predicted mean/cov; shift sigma points to the corrected mean (funcH linear -> valid).
   m_vecX = mean;
   m_matP = Pyy + m_matQmod;
   ensurePD(m_matP);
   for (int32_t i = 0; i < SIGMA_DIM_A; ++i)
      m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A).col(i) = sigmaY.col(i) + shift;

   // History for the smoother (one entry per predict, matching the standard path).
   m_vecXPredHist.push_back(m_vecX);
   m_matPPredHist.push_back(m_matP);
   // Smoother cross-corr C_i = Cov(x_filt_{i-1}, x_pred_i) = Pxx*F^T = Pxy
   // (NOT Pxy^T — Cov(a, F a) = E[a a^T] F^T = Pxx F^T). Matches the standard
   // path's calculateCrossCorrelation(input, filt, output, pred) = Cov(in,out).
   m_matCPredHist.push_back(Pxy);
}

void TrackFitterUKF::correctUKF(const ROOT::Math::XYZPoint &z)
{
   Vector<TF_DIM_Z> zVec; // Initialize the measurement vector
   zVec[0] = z.X();
   zVec[1] = z.Y();
   zVec[2] = z.Z();
   auto callback = [this](const kf::Vector<TF_DIM_X> &x_) { return funcH(x_); };
   TrackFitterUKFBase::correctUKF(callback, zVec);

   // logEigen("State PCorr", m_matP, m_matPCorrHist.size());
}

void TrackFitterUKF::smoothUKF()
{

   // Smoothing is done by iterating backwards over the history of predicted states and covariances.
   // Here i = k+1
   m_vecXSmooth.resize(m_vecXPredHist.size());
   m_matPSmooth.resize(m_matPPredHist.size());
   m_vecXSmooth.back() = m_vecXFiltHist.back(); // The last smoothed state is the last corrected state
   m_matPSmooth.back() = m_matPFiltHist.back(); // The last smoothed covariance is the last corrected covariance
   for (size_t i = m_vecXPredHist.size() - 1; i > 0; --i) {
      LOG(debug) << "Smoothing step " << i << " of " << m_vecXPredHist.size() - 1;

      // Get the predicted state and covariance at step i
      const auto &xPred = m_vecXPredHist[i]; // m_{k+1}^-
      const auto &pPred = m_matPPredHist[i]; // P_{k+1}^-
      const auto &ccor = m_matCPredHist[i];  // C_{k+1}

      // Get the filtered state and covariance at step i-1
      const auto &xFilt = m_vecXFiltHist[i - 1]; // m_{k}
      const auto &pFilt = m_matPFiltHist[i - 1]; // P_{k}

      // Get the smoothed state and covariance at step i
      auto &xSmooth = m_vecXSmooth[i]; // m^s_{k+1}
      auto &pSmooth = m_matPSmooth[i]; // P^s_{k+1}

      auto llt = calculateCholesky(pPred);
      auto D = ccor * llt.solve(Matrix<TF_DIM_X, TF_DIM_X>::Identity());
      // auto D = ccor * pPred.inverse(); // D = C_{k+1} * (P_{k+1}^-)^{-1}

      // std::cout << "D matrix at step " << i << ":\n" << D << "\n";
      m_vecXSmooth[i - 1] = xFilt + D * (xSmooth - xPred); // m^s_{k} = m_{k} + D * (m^s_{k+1} - m_{k+1}^-)
      m_matPSmooth[i - 1] =
         pFilt + D * (pSmooth - pPred) * D.transpose(); // P^s_{k} = P_{k} + D * (P^s_{k+1} - P_{k+1}^-) * D^T
      logEigen("State P+", m_matPSmooth[i - 1], i - 1);
   }
}

} // namespace kf