///
/// Copyright 2022 Mohanad Youssef (Al-khwarizmi)
///
/// Use of this source code is governed by an GPL-3.0 - style
/// license that can be found in the LICENSE file or at
/// https://opensource.org/licenses/GPL-3.0
///
/// @author Mohanad Youssef <mohanad.magdy.hammad@gmail.com>
/// @file unscented_kalman_filter.h
///

#ifndef TRACKFITTERUKF_H
#define TRACKFITTERUKF_H

#include "AtPropagator.h"

#include "kalman_filter.h"
#include "kf_util.h"
#include "unscented_kalman_filter.h"

namespace kf {

/// @brief Class for fitting tracks using the Unscented Kalman Filter (UKF) algorithm.
/// @tparam DIM_X Dimension of the state vector.
/// @tparam DIM_Z Dimension of the measurement vector.
/// @tparam DIM_V Dimension of the process noise vector.
/// @tparam DIM_N Dimension of the measurement noise vector.
template <int32_t DIM_X = 6, int32_t DIM_Z = 3, int32_t DIM_V = 2, int32_t DIM_N = 3>
class TrackFitterUKFBase : public KalmanFilter<DIM_X, DIM_Z> {
public:
   // Augmented state vector is just the process noise and state vector. The measurement noise is not included as that
   // is independent of the propagation and measurement model and just adds linearly.
   static constexpr int32_t DIM_A{DIM_X + DIM_V};       ///< @brief Augmented state dimension
   static constexpr int32_t SIGMA_DIM_A{2 * DIM_A + 1}; ///< @brief Sigma points dimension for augmented state
   float32_t m_kappa{0};                                ///< @brief Kappa parameter for finding sigma points

   // Add variables to track the covariances of the process and measurement noise.
   Matrix<DIM_V, DIM_V> m_matQ; // @brief Process noise covariance matrix
   Matrix<DIM_N, DIM_N> m_matR; // @brief Measurement noise covariance matrix

   Matrix<DIM_A, SIGMA_DIM_A> m_matSigmaXa{Matrix<DIM_A, SIGMA_DIM_A>::Zero()}; ///< @brief Sigma points matrix

   TrackFitterUKFBase()
      : KalmanFilter<DIM_X, DIM_Z>(), m_kappa(3 - DIM_A), m_matQ(Matrix<DIM_V, DIM_V>::Zero()),
        m_matR(Matrix<DIM_N, DIM_N>::Zero())
   {
      // 1. calculate weights
      updateWeights();
   }

   ~TrackFitterUKFBase() {}

   void setKappa(float32_t kappa)
   {
      m_kappa = kappa; // Set the kappa parameter for sigma point calculation
      updateWeights(); // Update the weights based on the new kappa value
   }

   // This code uses two different conventions for managing noise.
   // The state vector noise is set in the updateAugmentedStateAndCovariance() method, while
   // the noise vectors for the process and measurement models are set in the setCovarianceQ() and
   // setCovarianceR() methods. This is an odd choice. We will be moving everything into a common
   // structure where updateAugmentedStateAndCovariance() handle all covariance updates that are actually
   // part of the augmented state vector.

   ///
   /// @brief adding process noise covariance Q to the augmented state covariance
   /// matPa in the middle element of the diagonal.
   ///
   void setCovarianceQ(const Matrix<DIM_V, DIM_V> &matQ)
   {
      m_matQ = matQ; // Store the process noise covariance matrix
   }

   ///
   /// @brief set the measurement noise covariance R to be used in the update step
   ///
   void setCovarianceR(const Matrix<DIM_N, DIM_N> &matR)
   {
      m_matR = matR; // Store the measurement noise covariance matrix
   }

   /// Add state vector (m_vecX) to the augment state vector (m_vecXa) and also
   /// add covariance matrix (m_matP) to the augment covariance (m_matPa).
   void updateAugWithState()
   {
      // Copy state vector to augmented state vector
      for (int32_t i{0}; i < DIM_X; ++i) {
         m_vecXa[i] = m_vecX[i];
      }

      // Copy state covariance matrix to augmented covariance matrix
      for (int32_t i{0}; i < DIM_X; ++i) {
         for (int32_t j{0}; j < DIM_X; ++j) {
            m_matPa(i, j) = m_matP(i, j);
         }
      }
   }

   virtual std::array<float32_t, DIM_V> calculateProcessNoiseMean()
   {
      // Calculate the expectation value of the process noise using the current value of the state vector m_vecX
      std::array<float32_t, DIM_V> processNoiseMean{0};

      // TODO: Set the mean energy loss based on the momentum and particle type. Probably best to track stopping power?
      return processNoiseMean;
   }

   virtual Matrix<DIM_V, DIM_V> calculateProcessNoiseCovariance()
   {
      // Calculate the process noise covariance matrix
      Matrix<DIM_V, DIM_V> matQ{Matrix<DIM_V, DIM_V>::Zero()};
      matQ = m_matQ; // Use the stored process noise covariance matrix

      // TODO: Set the process noise covariance for angular straggle and energy loss.
      return matQ;
   }

   void updateAugWithProcessNoise()
   {
      auto processNoiseMean = calculateProcessNoiseMean();
      m_matQ = calculateProcessNoiseCovariance();

      // Add the mean process noise to the augmented state vector
      for (int32_t i{0}; i < DIM_V; ++i) {
         m_vecXa[DIM_X + i] = processNoiseMean[i];
      }

      // Add process noise covariance to the augmented covariance matrix
      const int32_t S_IDX{DIM_X};
      const int32_t L_IDX{S_IDX + DIM_V};

      for (int32_t i{S_IDX}; i < L_IDX; ++i) {
         for (int32_t j{S_IDX}; j < L_IDX; ++j) {
            m_matPa(i, j) = m_matQ(i - S_IDX, j - S_IDX);
         }
      }
   }

   ///
   /// @brief update the augmented state vector and covariance matrix
   /// This functions fully updates the augmented state vector (m_vecXa) and covariance matrix (m_matPa)
   /// by setting both the state vector and process noise components.
   ///
   void updateAugmentedStateAndCovariance()
   {
      updateAugWithState();
      updateAugWithProcessNoise();
   }

   ///
   /// @brief state prediction step of the unscented Kalman filter (UKF).
   /// @param predictionModelFunc callback to the prediction/process model
   /// function
   ///
   template <typename PredictionModelCallback>
   void predictUKF(PredictionModelCallback predictionModelFunc, const Vector<DIM_Z> &vecZ)
   {
      updateAugmentedStateAndCovariance();

      // Calculate the sigma points for the augmented state vector and save in a matrix where each column is a sigma
      // point.
      m_matSigmaXa = calculateSigmaPoints(m_vecXa, m_matPa);

      // Pull out the sigma points for the state vector and process noise in two different matrices.
      Matrix<DIM_X, SIGMA_DIM_A> sigmaXx{m_matSigmaXa.block(0, 0, DIM_X, SIGMA_DIM_A)}; // Sigma points for state vector
      Matrix<DIM_V, SIGMA_DIM_A> sigmaXv{
         m_matSigmaXa.block(DIM_X, 0, DIM_V, SIGMA_DIM_A)}; // Sigma points for process noise

      // Get each sigma point, apply the prediction model function, and store the results
      // back into the sigmaXx matrix (safe since each sigma point is independent).
      for (int32_t i{0}; i < SIGMA_DIM_A; ++i) {
         const Vector<DIM_X> sigmaXxi{util::getColumnAt<DIM_X, SIGMA_DIM_A>(i, sigmaXx)};
         const Vector<DIM_V> sigmaXvi{util::getColumnAt<DIM_V, SIGMA_DIM_A>(i, sigmaXv)};

         const Vector<DIM_X> Yi{predictionModelFunc(sigmaXxi, sigmaXvi, vecZ)}; // y = f(x)

         // Copy the predicted state vector back into the sigmaXx matrix.
         util::copyToColumn<DIM_X, SIGMA_DIM_A>(i, sigmaXx, Yi);
         // Copy the predicted state vector back into the augmented state vector (for use in future functions).
         util::copyToColumn<DIM_A, SIGMA_DIM_A, DIM_X>(i, m_matSigmaXa, Yi);
      }

      // Calculate the weighted mean and covariance of the sigma points for the state vector.
      // This will be the new state vector and covariance matrix.
      calculateWeightedMeanAndCovariance<DIM_X>(sigmaXx, m_vecX, m_matP);
   }

   ///
   /// @brief measurement correction step of the unscented Kalman filter (UKF).
   /// @param measurementModelFunc callback to the measurement model function
   /// @param vecZ actual measurement vector.
   ///
   template <typename MeasurementModelCallback>
   void correctUKF(MeasurementModelCallback measurementModelFunc, const Vector<DIM_Z> &vecZ)
   {
      // The state vector used here is an unaugmented state vector (m_vecX) and the covariance matrix is
      // an unaugmented covariance matrix (m_matP). This is because we are assuming the measurement noise
      // is independent of the state vector and process noise, so we can just use the unaugmented state vector
      // and covariance matrix for the measurement correction step, then add the measurement noise covariance.

      // Pull out the sigma points for the state vector after prediction.
      Matrix<DIM_X, SIGMA_DIM_A> sigmaXx{m_matSigmaXa.block(0, 0, DIM_X, SIGMA_DIM_A)}; // Sigma points for state vector

      // Get each sigma point, apply the prediction model function, and store the results
      // in the sigmaZ matrix.
      Matrix<DIM_Z, SIGMA_DIM_A> sigmaZ;
      for (int32_t i{0}; i < SIGMA_DIM_A; ++i) {
         const Vector<DIM_X> sigmaXxi{util::getColumnAt<DIM_X, SIGMA_DIM_A>(i, sigmaXx)};

         const Vector<DIM_Z> Zi{measurementModelFunc(sigmaXxi)}; // z = h(x)

         util::copyToColumn<DIM_Z, SIGMA_DIM_A>(i, sigmaZ, Zi);
      }

      // calculate the mean measurement vector and covariance matrix
      // from the sigma points.
      Vector<DIM_Z> vecZhat;
      Matrix<DIM_Z, DIM_Z> matPzz;
      calculateWeightedMeanAndCovariance<DIM_Z>(sigmaZ, vecZhat, matPzz);

      // Add in the measurement noise covariance matrix to the measurement covariance matrix.
      matPzz += m_matR; // Add measurement noise covariance

      // TODO: calculate cross correlation
      const Matrix<DIM_X, DIM_Z> matPxz{calculateCrossCorrelation(sigmaXx, m_vecX, sigmaZ, vecZhat)};

      // kalman gain
      const Matrix<DIM_X, DIM_Z> matK{matPxz * matPzz.inverse()};

      m_vecX += matK * (vecZ - vecZhat);
      m_matP -= matK * matPzz * matK.transpose();
   }

private:
   using KalmanFilter<DIM_X, DIM_Z>::m_vecX; // from Base KalmanFilter class
   using KalmanFilter<DIM_X, DIM_Z>::m_matP; // from Base KalmanFilter class

   float32_t m_weight0; /// @brief unscented transform weight 0 for mean
   float32_t m_weighti; /// @brief unscented transform weight i for none mean samples

   Vector<DIM_A> m_vecXa{Vector<DIM_A>::Zero()};               /// @brief augmented state vector (incl. process
                                                               /// and measurement noise means)
   Matrix<DIM_A, DIM_A> m_matPa{Matrix<DIM_A, DIM_A>::Zero()}; /// @brief augmented state covariance (incl.
                                                               /// process and measurement noise covariances)

   ///
   /// @brief algorithm to calculate the weights used to draw the sigma points
   ///
   void updateWeights()
   {
      static_assert(DIM_A > 0, "DIM_A is Zero which leads to numerical issue.");

      const float32_t denoTerm{m_kappa + static_cast<float32_t>(DIM_A)};

      m_weight0 = m_kappa / denoTerm;
      m_weighti = 0.5F / denoTerm;
   }

   ///
   /// @brief algorithm to calculate the deterministic sigma points for
   /// the unscented transformation
   ///
   /// @param vecX mean of the normally distributed state
   /// @param matPxx covariance of the normally distributed state
   /// @param STATE_DIM dimension of the vector used to calculate the sigma points
   /// @param SIGMA_DIM number of sigma points required (default is 2 * STATE_DIM + 1)
   /// @return matrix of sigma points where each column is a sigma point
   ///
   template <int32_t STATE_DIM, int32_t SIGMA_DIM = 2 * STATE_DIM + 1>
   Matrix<STATE_DIM, SIGMA_DIM>
   calculateSigmaPoints(const Vector<STATE_DIM> &vecXa, const Matrix<STATE_DIM, STATE_DIM> &matPa)
   {
      const float32_t scalarMultiplier{std::sqrt(STATE_DIM + m_kappa)}; // sqrt(n + \kappa)

      // cholesky factorization to get matrix Pxx square-root
      Eigen::LLT<Matrix<STATE_DIM, STATE_DIM>> lltOfPa(matPa);
      if (lltOfPa.info() != Eigen::Success) {
         throw std::runtime_error("Cholesky decomposition failed, matrix is not positive definite.");
      }
      Matrix<STATE_DIM, STATE_DIM> matSa{lltOfPa.matrixL()}; // sqrt(P_{a})

      matSa *= scalarMultiplier; // sqrt( (n + \kappa) * P_{a} )

      Matrix<STATE_DIM, SIGMA_DIM> sigmaXa;

      // X_0 = \bar{xa}
      util::copyToColumn<STATE_DIM, SIGMA_DIM>(0, sigmaXa, vecXa);

      for (int32_t i{0}; i < STATE_DIM; ++i) {
         const int32_t IDX_1{i + 1};
         const int32_t IDX_2{i + STATE_DIM + 1};

         util::copyToColumn<STATE_DIM, SIGMA_DIM>(IDX_1, sigmaXa, vecXa);
         util::copyToColumn<STATE_DIM, SIGMA_DIM>(IDX_2, sigmaXa, vecXa);

         const Vector<STATE_DIM> vecShiftTerm{util::getColumnAt<STATE_DIM, STATE_DIM>(i, matSa)};

         util::addColumnFrom<STATE_DIM, SIGMA_DIM>(IDX_1, sigmaXa, vecShiftTerm); // X_i^a     = \bar{xa} + sqrt( (n^a +
                                                                                  // \kappa) * P^{a} )
         util::subColumnFrom<STATE_DIM, SIGMA_DIM>(IDX_2, sigmaXa, vecShiftTerm); // X_{i+n}^a = \bar{xa} - sqrt( (n^a +
                                                                                  // \kappa) * P^{a} )
      }

      return sigmaXa;
   }

   ///
   /// @brief calculate the weighted mean and covariance given a set of sigma
   /// points
   /// @param[in] sigmaX matrix of (probably posterior) sigma points where each column contain single
   /// sigma point.
   /// @param[out] vecX output weighted mean of the sigma points
   /// @param[out] matPxx output weighted covariance of the sigma points
   ///
   template <int32_t STATE_DIM, int32_t SIGMA_DIM>
   void calculateWeightedMeanAndCovariance(const Matrix<STATE_DIM, SIGMA_DIM> &sigmaX, Vector<STATE_DIM> &vecX,
                                           Matrix<STATE_DIM, STATE_DIM> &matPxx)
   {
      // 1. calculate mean of the sigma points
      vecX = m_weight0 * util::getColumnAt<STATE_DIM, SIGMA_DIM>(0, sigmaX);
      for (int32_t i{1}; i < SIGMA_DIM; ++i) {
         vecX += m_weighti * util::getColumnAt<STATE_DIM, SIGMA_DIM>(i, sigmaX); // y += W[0, i] Y[:, i]
      }

      // 2. calculate covariance: P_{yy} = \sum_{i_0}^{2n} W[0, i] (Y[:, i] -
      // \bar{y}) (Y[:, i] - \bar{y})^T
      Vector<STATE_DIM> devXi{util::getColumnAt<STATE_DIM, SIGMA_DIM>(0, sigmaX) - vecX}; // Y[:, 0] - \bar{ y }

      matPxx = m_weight0 * devXi * devXi.transpose(); // P_0 = W[0, 0] (Y[:, 0] - \bar{y}) (Y[:, 0] -
                                                      // \bar{y})^T

      for (int32_t i{1}; i < SIGMA_DIM; ++i) {
         devXi = util::getColumnAt<STATE_DIM, SIGMA_DIM>(i, sigmaX) - vecX; // Y[:, i] - \bar{y}

         const Matrix<STATE_DIM, STATE_DIM> Pi{m_weighti * devXi * devXi.transpose()}; // P_i = W[0, i] (Y[:, i] -
                                                                                       // \bar{y}) (Y[:, i] - \bar{y})^T

         matPxx += Pi; // y += W[0, i] (Y[:, i] - \bar{y}) (Y[:, i] - \bar{y})^T
      }
   }

   ///
   /// @brief calculate the cross-correlation given two sets sigma points X and Y
   /// and their means x and y
   /// @param sigmaX first matrix of sigma points where each column contain
   /// single sigma point
   /// @param vecX mean of the first set of sigma points
   /// @param sigmaY second matrix of sigma points where each column contain
   /// single sigma point
   /// @param vecY mean of the second set of sigma points
   /// @return matPxy, the cross-correlation matrix
   ///
   template <int32_t SIGMA_DIM>
   Matrix<DIM_X, DIM_Z> calculateCrossCorrelation(const Matrix<DIM_X, SIGMA_DIM> &sigmaX, const Vector<DIM_X> &vecX,
                                                  const Matrix<DIM_Z, SIGMA_DIM> &sigmaY, const Vector<DIM_Z> &vecY)
   {
      Vector<DIM_X> devXi{util::getColumnAt<DIM_X, SIGMA_DIM>(0, sigmaX) - vecX}; // X[:, 0] - \bar{ x }
      Vector<DIM_Z> devYi{util::getColumnAt<DIM_Z, SIGMA_DIM>(0, sigmaY) - vecY}; // Y[:, 0] - \bar{ y }

      // P_0 = W[0, 0] (X[:, 0] - \bar{x}) (Y[:, 0] - \bar{y})^T
      Matrix<DIM_X, DIM_Z> matPxy{m_weight0 * (devXi * devYi.transpose())};

      for (int32_t i{1}; i < SIGMA_DIM; ++i) {
         devXi = util::getColumnAt<DIM_X, SIGMA_DIM>(i, sigmaX) - vecX; // X[:, i] - \bar{x}
         devYi = util::getColumnAt<DIM_Z, SIGMA_DIM>(i, sigmaY) - vecY; // Y[:, i] - \bar{y}

         matPxy += m_weighti * (devXi * devYi.transpose()); // y += W[0, i] (Y[:, i] -
                                                            // \bar{y}) (Y[:, i] - \bar{y})^T
      }

      return matPxy;
   }
};

// Define template dimension variables for clarity and reuse

class TrackFitterUKF : public TrackFitterUKFBase<6, 1, 3, 3> {
protected:
   static constexpr int32_t TRACKFITTER_DIM_X = 6;
   static constexpr int32_t TRACKFITTER_DIM_Z = 1;
   static constexpr int32_t TRACKFITTER_DIM_V = 3;
   static constexpr int32_t TRACKFITTER_DIM_N = 3;

   AtTools::AtPropagator fPropagator;          ///< @brief Propagator for the track fitter
   AtTools::AtPropagator::StepState fMeanStep; /// Holds the step information for POCA propagation of mean state

public:
   /**
    * @brief Constructor for the TrackFitterUKF class.
    * @param propagator The propagator to be used for the track fitting, must be passed as an rvalue reference.
    *
    * Example usage:
    * ```
    * AtTools::AtPropagator propagator;
    * kf::TrackFitterUKF trackFitterUKF(std::move(propagator));
    * ```
    */
   TrackFitterUKF(AtTools::AtPropagator &&propagator) : TrackFitterUKFBase(), fPropagator(std::move(propagator)) {}

   template <typename PredictionModelCallback>
   void predictUKF(PredictionModelCallback predictionModelFunc, const Vector<TRACKFITTER_DIM_Z> &vecZ)
   {
      // First we need to propagate the mean state vector to the next measurement point.
      AtTools::AtRK4Stepper stepper; // Use RK4 stepper for propagation
      fPropagator.PropagateToMeasurementSurface(AtTools::AtMeasurementPoint(vecZ), stepper);
      fMeanStep = fPropagator.GetState(); // Get the mean step information from the propagator
   }

protected:
   std::array<float32_t, TRACKFITTER_DIM_V> calculateProcessNoiseMean() override
   {
      // The process noise is the scaling factor for dedx. By definition the mean should be 1
      std::array<float32_t, TRACKFITTER_DIM_V> processNoiseMean{1};
      return processNoiseMean;
   }

   Matrix<TRACKFITTER_DIM_V, TRACKFITTER_DIM_V> calculateProcessNoiseCovariance() override
   {
      assert(TRACKFITTER_DIM_V == 1 && "Process noise covariance is only implemented for DIM_V = 1");
      // Calculate the process noise covariance matrix
      Matrix<TRACKFITTER_DIM_V, TRACKFITTER_DIM_V> matQ{Matrix<TRACKFITTER_DIM_V, TRACKFITTER_DIM_V>::Zero()};

      // We need to know what the energy of the particle before/after transport.
      double eIn = 0;
      double eOut = 0;

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
};
} // namespace kf

#endif // TRACKFITTERUKF_H
