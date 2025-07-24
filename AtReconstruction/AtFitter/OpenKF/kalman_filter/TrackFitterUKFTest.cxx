#include "kalman_filter/TrackFitterUKF.h"

#include "AtELossTable.h"
#include "AtKinematics.h"
#include "AtPropagator.h"

#include <FairLogger.h>

#include "TMatrixD.h"

#include "Math/Vector3D.h"
#include "gtest/gtest.h"

class TrackFitterUKFExampleTest : public testing::Test {
public:
   virtual void SetUp() override {}
   virtual void TearDown() override {}

   static constexpr float FLOAT_EPSILON{0.001F};

   static constexpr size_t DIM_X{4};
   static constexpr size_t DIM_V{4};
   static constexpr size_t DIM_Z{2};
   static constexpr size_t DIM_N{2};

   kf::TrackFitterUKFBase<DIM_X, DIM_Z, DIM_V, DIM_N> m_ukf;

   /// @brief to propagate the state vector using the process model
   /// @param x state vector
   /// @param v process noise vector
   /// @return propagated (unaugmented) state vector
   static kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> &x, const kf::Vector<DIM_V> &v, const kf::Vector<DIM_Z> &z)
   {
      kf::Vector<DIM_X> y;
      y[0] = x[0] + x[2] + v[0];
      y[1] = x[1] + x[3] + v[1];
      y[2] = x[2] + v[2];
      y[3] = x[3] + v[3];
      return y;
   }

   /// @brief to apply the measurement model to the state vector
   /// @param x the state vector of the system
   /// @return the measurement vector
   static kf::Vector<DIM_Z> funcH(const kf::Vector<DIM_X> &x)
   {
      kf::Vector<DIM_Z> y;

      kf::float32_t px{x[0]};
      kf::float32_t py{x[1]};

      y[0] = std::sqrt((px * px) + (py * py));
      y[1] = std::atan(py / (px + std::numeric_limits<kf::float32_t>::epsilon()));
      return y;
   }
};

TEST_F(TrackFitterUKFExampleTest, Prediction)
{
   kf::Vector<DIM_X> x;
   x << 2.0F, 1.0F, 0.0F, 0.0F;

   kf::Matrix<DIM_X, DIM_X> P;
   P << 0.01F, 0.0F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F;

   kf::Matrix<DIM_V, DIM_V> Q;
   Q << 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.1F, 0.0F, 0.0F, 0.0F, 0.0F, 0.1F;

   kf::Matrix<DIM_N, DIM_N> R;
   R << 0.01F, 0.0F, 0.0F, 0.01F;

   kf::Vector<DIM_Z> z;
   z << 2.5F, 0.05F;

   m_ukf.vecX() = x;
   m_ukf.matP() = P;

   m_ukf.setCovarianceQ(Q);
   m_ukf.setCovarianceR(R);

   m_ukf.predictUKF(funcF, z);

   // Expectation from the python results:
   // =====================================
   // x =
   //     [2.0 1.0 0.0 0.0]
   // P =
   //     [[0.11  0.00  0.05  0.00]
   //      [0.00  0.11  0.00  0.05]
   //      [0.05  0.00  0.15  0.00]
   //      [0.00  0.05  0.00  0.15]]

   ASSERT_NEAR(m_ukf.vecX()[0], 2.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[1], 1.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[2], 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[3], 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(0, 0), 0.11F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 1), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 2), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(1, 0), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 1), 0.11F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 3), 0.05F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(2, 0), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 1), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 2), 0.15F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(3, 0), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 1), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 3), 0.15F, FLOAT_EPSILON);
}

TEST_F(TrackFitterUKFExampleTest, PredictionAndCorrection)
{
   kf::Vector<DIM_X> x;
   x << 2.0F, 1.0F, 0.0F, 0.0F;

   kf::Matrix<DIM_X, DIM_X> P;
   P << 0.01F, 0.0F, 0.0F, 0.0F, 0.0F, 0.01F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F;

   kf::Matrix<DIM_V, DIM_V> Q;
   Q << 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.05F, 0.0F, 0.0F, 0.0F, 0.0F, 0.1F, 0.0F, 0.0F, 0.0F, 0.0F, 0.1F;

   kf::Matrix<DIM_N, DIM_N> R;
   R << 0.01F, 0.0F, 0.0F, 0.01F;

   kf::Vector<DIM_Z> z;
   z << 2.5F, 0.05F;

   m_ukf.vecX() = x;
   m_ukf.matP() = P;

   m_ukf.setCovarianceQ(Q);
   m_ukf.setCovarianceR(R);

   m_ukf.predictUKF(funcF, z);

   // Expectation from the python results:
   // =====================================
   // x =
   //     [2.0 1.0 0.0 0.0]
   // P =
   //     [[0.11  0.00  0.05  0.00]
   //      [0.00  0.11  0.00  0.05]
   //      [0.05  0.00  0.15  0.00]
   //      [0.00  0.05  0.00  0.15]]

   ASSERT_NEAR(m_ukf.vecX()[0], 2.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[1], 1.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[2], 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[3], 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(0, 0), 0.11F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 1), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 2), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(1, 0), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 1), 0.11F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 3), 0.05F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(2, 0), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 1), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 2), 0.15F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(3, 0), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 1), 0.05F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 3), 0.15F, FLOAT_EPSILON);

   m_ukf.correctUKF(funcH, z);

   // Expectations from the python results:
   // ======================================
   // x =
   //     [ 2.4758845   0.53327217  0.21649734 -0.21214576]
   // P =
   // [[ 0.01433114 -0.01026142  0.00651178 -0.00465059]
   //  [-0.01026142  0.0295458  -0.0046378   0.01344241]
   //  [ 0.00651178 -0.0046378   0.13023154 -0.00210188]
   //  [-0.00465059  0.01344241 -0.00210188  0.1333886 ]]

   ASSERT_NEAR(m_ukf.vecX()[0], 2.4758845F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[1], 0.53327217F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[2], 0.21649734F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[3], -0.21214576F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(0, 0), 0.01433114F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 1), -0.01026142F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 2), 0.00651178F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 3), -0.00465059F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(1, 0), -0.01026142F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 1), 0.0295458F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 2), -0.0046378F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 3), 0.01344241F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(2, 0), 0.00651178F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 1), -0.0046378F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 2), 0.13023154F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 3), -0.00210188F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(3, 0), -0.00465059F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 1), 0.01344241F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 2), -0.00210188F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 3), 0.1333886F, FLOAT_EPSILON);
}

namespace {

std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr) {
      return "../../resources/energy_loss/HinH.txt"; // Default path assuming cwd is build/AtTools
   }
   return std::string(env) + "/resources/energy_loss/HinH.txt"; // Use environment variable
}
} // namespace

namespace kf {
class TrackFitterUKFTest : public TrackFitterUKF {
public:
   Vector<DIM_A> vecXa() { return m_vecXa; }
   Matrix<DIM_A, DIM_A> matPa() { return m_matPa; }
   Matrix<DIM_A, SIGMA_DIM_A> matSigmaXa() { return m_matSigmaXa; }
   AtTools::AtPropagator &getPropagator() { return fPropagator; }
   AtTools::AtPropagator::StepState &getState() { return fMeanStep; }

   void calcSigmaPoints() { m_matSigmaXa = calculateSigmaPoints(m_vecXa, m_matPa); }
   void calcMeanAndCovFromSigma(const Matrix<DIM_A, SIGMA_DIM_A> &sigmaXa, Vector<DIM_A> &meanXa,
                                Matrix<DIM_A, DIM_A> &covXa)
   {
      calculateWeightedMeanAndCovariance(sigmaXa, meanXa, covXa);
   }

   TrackFitterUKFTest(AtTools::AtPropagator &&propagator, std::unique_ptr<AtTools::AtStepper> &&stepper)
      : TrackFitterUKF(std::move(propagator), std::move(stepper))
   {
   }
};
} // namespace kf
class TrackFitterUKFFixture : public testing::Test {
public:
   using XYZVector = ROOT::Math::XYZVector;

   virtual void TearDown() override {}

   static constexpr float FLOAT_EPSILON{0.001F};

   static constexpr size_t DIM_X{6}; // Set state vector to be (x,y,z,p,theta,phi)
   static constexpr size_t DIM_V{1}; // Noise vector is energy straggling
   static constexpr size_t DIM_Z{3}; // Measurement vector is (x,y,z)
   static constexpr size_t DIM_N{3}; // Measurement noise vector

   XYZVector fBField{0, 0, 2.85};           // B-field in tesla
   XYZVector fEField{0, 0, 0};              // E-field in V/m
   static constexpr double fC = 299792458;  // m/s
   const double mass_p = 938.272;           // Mass of proton in MeV/c^2
   const double charge_p = 1.602176634e-19; // Charge of protonble

   std::unique_ptr<kf::TrackFitterUKF> m_ukf{nullptr};

   void SetUp() override
   {
      auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
      elossModel->LoadSrimTable(getEnergyPath()); // Use the function to get the path
      AtTools::AtPropagator propagator(charge_p, mass_p, std::move(elossModel));
      propagator.SetBField(fBField);
      propagator.SetEField(fEField);
      auto stepper = std::make_unique<AtTools::AtRK4Stepper>(); // Use RK4 stepper for propagation
      m_ukf = std::make_unique<kf::TrackFitterUKFTest>(std::move(propagator), std::move(stepper));
   }

   /// @brief to propagate the state vector using the process model
   /// @param x state vector
   /// @param v process noise vector
   /// @param vecZ The next measurement point used to stop the propagation.
   /// @return propagated (unaugmented) state vector
   static kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> &x, const kf::Vector<DIM_V> &v, const kf::Vector<DIM_Z> &z)
   {
      return {};
   }
};

TEST_F(TrackFitterUKFFixture, TestInstantiation)
{
   assert(true);
}

TEST_F(TrackFitterUKFFixture, TestAugmentedStateUpdate)
{
   auto testPtr = dynamic_cast<kf::TrackFitterUKFTest *>(m_ukf.get());

   // Create a state change where a proton has lost 1 MeV of energy (to rest)
   auto p = AtTools::Kinematics::GetRelMomFromKE(1.0, mass_p);
   AtTools::AtPropagator::StepState fInitialState;
   fInitialState.fMass = mass_p;
   fInitialState.fQ = charge_p;

   // Current state
   fInitialState.fLastPos = {0, 0, 0};
   fInitialState.fLastMom = {p, 0, 0};
   // State at next measurement point
   fInitialState.fPos = {202.78, 0, 0}; // Pulled from SRIM table
   fInitialState.fMom = {0, 0, 0};
   double var = 0.01;
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 6; ++i) {
      cov(i, i) = var; // Set diagonal covariance to some small number
   }

   m_ukf->SetInitialState(fInitialState.fLastPos, fInitialState.fLastMom, cov);
   testPtr->getState() = fInitialState; // Set the mean state in the test class (rather than propagate)

   // Fill m_vecXa and m_matPa with the provided state
   m_ukf->updateAugmentedStateAndCovariance();

   assert(testPtr != nullptr);

   // Check the augmented state vector
   double eps = 1e-4;
   ASSERT_NEAR(testPtr->vecXa()[0], fInitialState.fLastPos.X(), eps);
   ASSERT_NEAR(testPtr->vecXa()[1], fInitialState.fLastPos.Y(), eps);
   ASSERT_NEAR(testPtr->vecXa()[2], fInitialState.fLastPos.Z(), eps);
   ASSERT_NEAR(testPtr->vecXa()[3], fInitialState.fLastMom.R(), eps);
   ASSERT_NEAR(testPtr->vecXa()[4], fInitialState.fPos.Theta(), eps);
   ASSERT_NEAR(testPtr->vecXa()[5], fInitialState.fPos.Phi(), eps);
   ASSERT_NEAR(testPtr->vecXa()[6], 1, eps); // Average energy straggling factor is 1

   // Check the augmented covariance matrix. Nothing set but energy straggling
   ASSERT_NEAR(testPtr->matPa()(0, 0), var, eps);
   ASSERT_NEAR(testPtr->matPa()(1, 1), var, eps);
   ASSERT_NEAR(testPtr->matPa()(2, 2), var, eps);
   ASSERT_NEAR(testPtr->matPa()(3, 3), var, eps);
   ASSERT_NEAR(testPtr->matPa()(4, 4), var, eps);
   ASSERT_NEAR(testPtr->matPa()(5, 5), var, eps);

   double sigma = 8.74 / 202.78; // Energy straggling factor (from range straggling table)
   ASSERT_NEAR(testPtr->matPa()(6, 6), sigma * sigma, eps);
}

TEST_F(TrackFitterUKFFixture, TestSigmaPoints)
{
   auto testPtr = dynamic_cast<kf::TrackFitterUKFTest *>(m_ukf.get());

   // Create a state change where a proton has lost 1 MeV of energy (to rest)
   auto p = AtTools::Kinematics::GetRelMomFromKE(1.0, mass_p);
   AtTools::AtPropagator::StepState fInitialState;
   fInitialState.fMass = mass_p;
   fInitialState.fQ = charge_p;

   // Current state
   fInitialState.fLastPos = {0, 0, 0};
   fInitialState.fLastMom = {p, 0, 0};
   // State at next measurement point
   fInitialState.fPos = {202.78, 0, 0}; // Pulled from SRIM table
   fInitialState.fMom = {0, 0, 0};
   double var = 0.01;
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 6; ++i) {
      cov(i, i) = var; // Set diagonal covariance to some small number
   }

   m_ukf->SetInitialState(fInitialState.fLastPos, fInitialState.fLastMom, cov);
   testPtr->getState() = fInitialState; // Set the mean state in the test class (rather than propagate)

   // Fill m_vecXa and m_matPa with the provided state
   m_ukf->updateAugmentedStateAndCovariance();

   testPtr->calcSigmaPoints();
   auto sigmaPoints = testPtr->matSigmaXa();

   kf::Matrix<kf::TrackFitterUKF::DIM_A, kf::TrackFitterUKF::DIM_A> covXa;
   kf::Vector<kf::TrackFitterUKF::DIM_A> meanXa;
   testPtr->calcMeanAndCovFromSigma(sigmaPoints, meanXa, covXa);

   // Verify the mean and cov match what we put in
   for (int i = 0; i < kf::TrackFitterUKF::DIM_A; ++i) {
      ASSERT_NEAR(meanXa[i], testPtr->vecXa()[i], FLOAT_EPSILON);
      for (int j = 0; j < kf::TrackFitterUKF::DIM_A; ++j) {
         ASSERT_NEAR(covXa(i, j), testPtr->matPa()(i, j), FLOAT_EPSILON);
      }
   }

   LOG(debug) << "Sigma points:\n" << sigmaPoints;
}

TEST_F(TrackFitterUKFFixture, TestPredictionStep)
{
   auto testPtr = dynamic_cast<kf::TrackFitterUKFTest *>(m_ukf.get());

   // Create a state change where a proton has lost 1 MeV of energy (to rest)
   auto p = AtTools::Kinematics::GetRelMomFromKE(1.0, mass_p);
   AtTools::AtPropagator::StepState fInitialState;
   fInitialState.fMass = mass_p;
   fInitialState.fQ = charge_p;

   // Current state
   fInitialState.fLastPos = {0, 0, 0};
   fInitialState.fLastMom = {p, 0, 0};
   // State at next measurement point
   fInitialState.fPos = {202.78, 0, 0}; // Pulled from SRIM table
   fInitialState.fMom = {0, 0, 0};
   double var = 0.01;
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 6; ++i) {
      cov(i, i) = var; // Set diagonal covariance to some small number
   }

   m_ukf->SetInitialState(fInitialState.fLastPos, fInitialState.fLastMom, cov);
   testPtr->getState() = fInitialState; // Set the mean state in the test class (rather than propagate)

   // Fill m_vecXa and m_matPa with the provided state
   m_ukf->updateAugmentedStateAndCovariance();
}