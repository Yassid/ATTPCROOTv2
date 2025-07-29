#include "TrackFitterUKF.h"

#include <FairLogger.h>

#include <Math/Plane3D.h>
#include <TMatrixD.h>

#include <AtPropagator.h>

namespace kf {
void TrackFitterUKF::Reset()
{
   // Reset the state vector and covariance matrix
   m_vecX.setZero();
   m_matP.setZero();
   m_vecXa.setZero();
   m_matPa.setZero();
   m_matQ.setZero();
   m_matR.setZero();
   m_matSigmaXa.setZero();

   // Clear the history vectors
   m_vecXPredHist.clear();
   m_matPPredHist.clear();
   m_matCPredHist.clear();
   m_vecXHist.clear();
   m_matPHist.clear();
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

   // Save the initial state in our history vectors
   m_vecXHist.push_back(m_vecX);
   m_matPHist.push_back(m_matP);
   m_vecXPredHist.push_back(m_vecX);
   m_matPPredHist.push_back(m_matP);
   m_matCPredHist.push_back(Matrix<TF_DIM_X, TF_DIM_X>::Zero()); // Cross-correlation is not defined for the first point

   // We need to calculate the sigma points for the initial state
   updateAugmentedStateAndCovariance();                   // Update the augmented state vector and covariance matrix
   m_matSigmaXa = calculateSigmaPoints(m_vecXa, m_matPa); // Calculate the sigma points for the initial state
   // Now we grab the sigma points only for the state.
   m_matSigmaXPred = m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A); // Extract the state sigma points
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
   return {m_vecXa[0], m_vecXa[1], m_vecXa[2], m_vecXa[3],
           m_vecXa[4], m_vecXa[5], m_vecXa[6]}; // Return the augmented state vector as an array
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
      // If energy straggling is disabled, we set the process noise to zero.
      matQ(0, 0) = 0.0;
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
      LOG(info) << "Calculating process noise for straggling between " << eIn << " MeV and " << eOut << " MeV over "
                << elossModel->GetRange(eIn, eOut) << " mm.";
      LOG(info) << "Process noise covariance for energy straggling: " << matQ(0, 0) << " (factor: " << factor
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
   XYZPoint fPos(x[0], x[1], x[2]);      // Position from state vector
   Polar3DVector fMom(x[3], x[4], x[5]); // Momentum from state vector

   fPropagator.SetState(fPos, XYZVector(fMom));

   fPropagator.fScalingFactor = v[0]; // Set the scaling factor for energy loss
   fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPlane(fMeasurementPlane), *fStepper);
   fPropagator.fScalingFactor = 1.0; // Reset the scaling factor after propagation

   auto fState = fPropagator.GetState(); // Get the propagated state
   Vector<TF_DIM_X> vecX{Vector<TF_DIM_X>::Zero()};
   vecX[0] = fState.fPos.X();     // X position
   vecX[1] = fState.fPos.Y();     // Y position
   vecX[2] = fState.fPos.Z();     // Z position
   vecX[3] = fState.fMom.R();     // Momentum magnitude
   vecX[4] = fState.fMom.Theta(); // Polar angle
   vecX[5] = fState.fMom.Phi();   // Azimuthal

   return vecX; // Return the propagated state vector
}

Vector<TrackFitterUKF::TF_DIM_Z> TrackFitterUKF::funcH(const Vector<TrackFitterUKF::TF_DIM_X> &x)
{
   // Measurement function to convert state vector to measurement vector
   using namespace ROOT::Math;
   Vector<TF_DIM_Z> vecZ;
   XYZPoint fPos(x[0], x[1], x[2]); // Position from state vector

   // Calculate the measurement vector based on the position and momentum
   vecZ[0] = fPos.X(); // X coordinate
   vecZ[1] = fPos.Y(); // Y coordinate
   vecZ[2] = fPos.Z(); // Z coordinate

   return vecZ; // Return the measurement vector
}

void TrackFitterUKF::predictUKF(const ROOT::Math::XYZPoint &z)
{
   using namespace ROOT::Math;

   // First we need to propagate the mean state vector to the next measurement point.
   XYZPoint startingPosition{m_vecX[0], m_vecX[1], m_vecX[2]};      // Get the starting position from the state vector
   Polar3DVector startingMomentum{m_vecX[3], m_vecX[4], m_vecX[5]}; // Get the starting momentum from the state vector

   LOG(info) << "Propagating reference state from position: " << startingPosition
             << " with momentum: " << XYZVector(startingMomentum);

   fPropagator.SetState(startingPosition, XYZVector(startingMomentum));
   fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPoint(z), *fStepper);
   fMeanStep = fPropagator.GetState();    // Get the mean step information from the propagator
   fMeanStep.fLastPos = startingPosition; // Store the last position
   fMeanStep.fLastMom = startingMomentum; // Store the last momentum

   LOG(info) << "Propagated to position: " << fMeanStep.fPos << " with momentum: " << fMeanStep.fMom;

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

   // Now we need to store the predicted state and covariance for smoothing later.
   m_vecXPredHist.push_back(m_vecX); // Store the predicted state vector
   m_matPPredHist.push_back(m_matP); // Store the predicted covariance matrix

   // Get the sigma points belonging to the predicted state
   Matrix<TF_DIM_X, SIGMA_DIM_A> sigmaXx{m_matSigmaXa.block(0, 0, TF_DIM_X, SIGMA_DIM_A)};

   // Calculate the cross-corelation between the filtered state at k and predicted state at k+1
   auto matCPred =
      calculateCrossCorrelation<TF_DIM_X>(m_matSigmaXPred, m_vecXHist.back(), sigmaXx, m_vecXPredHist.back());
   m_matCPredHist.push_back(matCPred); // Store the cross-correlation matrix
}

void TrackFitterUKF::correctUKF(const ROOT::Math::XYZPoint &z)
{
   Vector<TF_DIM_Z> zVec; // Initialize the measurement vector
   zVec[0] = z.X();
   zVec[1] = z.Y();
   zVec[2] = z.Z();
   auto callback = [this](const kf::Vector<TF_DIM_X> &x_) { return funcH(x_); };
   TrackFitterUKFBase::correctUKF(callback, zVec);

   // After correction we need to save the filtered state
   m_vecXHist.push_back(m_vecX); // Store the filtered state vector
   m_matPHist.push_back(m_matP); // Store the filtered covariance matrix
}

} // namespace kf