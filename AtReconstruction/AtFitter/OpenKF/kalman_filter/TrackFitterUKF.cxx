#include "TrackFitterUKF.h"

#include <FairLogger.h>

#include <Math/Plane3D.h>
#include <TMatrixD.h>

#include <AtPropagator.h>

namespace kf {

void TrackFitterUKF::SetInitialState(const ROOT::Math::XYZPoint &initialPosition,
                                     const ROOT::Math::XYZVector &initialMomentum, const TMatrixD &initialCovariance)
{
   m_vecX[0] = initialPosition.X();     // X position
   m_vecX[1] = initialPosition.Y();     // Y position
   m_vecX[2] = initialPosition.Z();     // Z position
   m_vecX[3] = initialMomentum.R();     // Momentum magnitude
   m_vecX[4] = initialMomentum.Theta(); // Polar angle
   m_vecX[5] = initialMomentum.Phi();   // Azimuthal angle

   // Copy elements from initialCovariance to m_matP
   for (int i = 0; i < m_matP.rows(); ++i) {
      for (int j = 0; j < m_matP.cols(); ++j) {
         m_matP(i, j) = initialCovariance(i, j);
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

   if (const auto *elossModel = fPropagator.GetELossModel()) {
      double dedx_straggle = elossModel->GetdEdxStraggling(eIn, eOut);
      double factor = dedx_straggle / elossModel->GetdEdx(eIn);
      matQ(0, 0) = factor * factor; // Variance for the dedx straggling.

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

} // namespace kf