#include "kalman_filter/AtTrackFitterUKF.h"

#include "gtest/gtest.h"

class AtTrackFitterUKFTest : public testing::Test {
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

TEST_F(AtTrackFitterUKFTest, test_UKF_Prediction)
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

TEST_F(AtTrackFitterUKFTest, test_UKF_PredictionAndCorrection)
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
   //     [ 2.554  0.356  0.252 -0.293]
   // P =
   //     [[ 0.01  -0.001  0.005 -0.    ]
   //      [-0.001  0.01 - 0.     0.005 ]
   //      [ 0.005 - 0.     0.129 - 0.  ]
   //      [-0.     0.005 - 0.     0.129]]

   ASSERT_NEAR(m_ukf.vecX()[0], 2.554F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[1], 0.356F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[2], 0.252F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.vecX()[3], -0.293F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(0, 0), 0.01F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 1), -0.001F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 2), 0.005F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(0, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(1, 0), -0.001F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 1), 0.01F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(1, 3), 0.005F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(2, 0), 0.005F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 1), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 2), 0.129F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(2, 3), 0.0F, FLOAT_EPSILON);

   ASSERT_NEAR(m_ukf.matP()(3, 0), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 1), 0.005F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 2), 0.0F, FLOAT_EPSILON);
   ASSERT_NEAR(m_ukf.matP()(3, 3), 0.129F, FLOAT_EPSILON);
}
