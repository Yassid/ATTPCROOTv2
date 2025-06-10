#include "kalman_filter/TrackFitterUKF.h"

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

   kf::TrackFitterUKF<DIM_X, DIM_Z, DIM_V, DIM_N> m_ukf;

   /// @brief to propagate the state vector using the process model
   /// @param x state vector
   /// @param v process noise vector
   /// @return propagated (unaugmented) state vector
   static kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> &x, const kf::Vector<DIM_V> &v)
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

   m_ukf.predictUKF(funcF);

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

   m_ukf.predictUKF(funcF);

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

class TrackFitterUKFPhysicsTest : public testing::Test {
public:
   virtual void SetUp() override {}
   virtual void TearDown() override {}

   static constexpr float FLOAT_EPSILON{0.001F};

   static constexpr size_t DIM_X{6};
   static constexpr size_t DIM_V{2};
   static constexpr size_t DIM_Z{3};
   static constexpr size_t DIM_N{3};

   kf::TrackFitterUKF<DIM_X, DIM_Z, DIM_V, DIM_N> m_ukf;

   /// @brief to propagate the state vector using the process model
   /// @param x state vector
   /// @param v process noise vector
   /// @return propagated (unaugmented) state vector
   static kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> &x, const kf::Vector<DIM_V> &v)
   {
      //TODO: This needs to be filled with an RK4 solver for the physics model
      kf::Vector<DIM_X> y{x};
      
      // For now, we just return the state vector as is
      return y;
   }

   /// @brief to apply the measurement model to the state vector
   /// @param x the state vector of the system
   /// @return the measurement vector
   static kf::Vector<DIM_Z> funcH(const kf::Vector<DIM_X> &x)
   {
      kf::Vector<DIM_Z> y;
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2];

      return y;
   }
};

TEST_F(TrackFitterUKFPhysicsTest, PhysicsPrediction)
{
   // TODO: This needs to be filled with a proper physics test. 
   kf::Vector<DIM_X> x; // Initial state vector

   kf::Matrix<DIM_X, DIM_X> P; // Initial state vector covariance matrix

   // Note: process noise is defined in the UKF class, so we don't need to set it here.

   kf::Vector<DIM_Z> z; // Measurement vector to be used in the correction step
   z << std::sqrt(5), std::atan2(1.f ,2.f);

   kf::Matrix<DIM_N, DIM_N> R; // Covariance matrix for the measurement noise

   m_ukf.vecX() = x;
   m_ukf.matP() = P;

   m_ukf.setCovarianceR(R);
   m_ukf.predictUKF(funcF);


}
