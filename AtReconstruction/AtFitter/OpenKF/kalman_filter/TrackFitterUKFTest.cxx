#include "kalman_filter/TrackFitterUKF.h"

#include "AtELossTable.h"
#include "AtKinematics.h"

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

   kf::TrackFitterUKF<DIM_X, DIM_Z, DIM_V, DIM_N> m_ukf;

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

class TrackFitterUKFPhysicsTest : public testing::Test {
public:
   using XYZVector = ROOT::Math::XYZVector;

   virtual void SetUp() override {}
   virtual void TearDown() override {}

   static constexpr float FLOAT_EPSILON{0.001F};

   static constexpr size_t DIM_X{6};
   static constexpr size_t DIM_V{2};
   static constexpr size_t DIM_Z{3};
   static constexpr size_t DIM_N{3};

   static XYZVector fBField;               // B-field in tesla
   static XYZVector fEField;               // E-field in V/m
   static constexpr double fC = 299792458; // m/s

   kf::TrackFitterUKF<DIM_X, DIM_Z, DIM_V, DIM_N> m_ukf;

   /**
    * @param pos Position of particle in mm
    * @param mom Momentum of particle in MeV/c
    * @param charge charge of the particle in Coulombs
    * @param mass mass of particle in MeV/c^2
    * @param dedx Stopping power in MeV/mm
    *
    * @returns Force in N
    */
   static XYZVector Force(XYZVector pos, XYZVector mom, double charge, double mass, double dedx)
   {

      // auto fourMom = AtTools::Kinematics::Get4Vector(mom, mass);
      // auto v = mom / fourMom.E() * c; // m/s
      auto v = GetVel(mom, mass);

      auto F_lorentz = charge * (fEField + v.Cross(fBField));
      // std::cout << "F_lorentz: " << F_lorentz << std::endl;
      auto dedx_si = dedx * 1.60218e-10; // de_dx in SI units (J/m)

      auto drag = -dedx_si * mom.Unit();
      // std::cout << "drag: " << drag << " mom " << mom << " dedx " << dedx_si << std::endl;

      return F_lorentz + drag;
   }

   static XYZVector GetVel(XYZVector mom, double mass)
   {
      auto fourMom = AtTools::Kinematics::Get4Vector(mom, mass);
      const double c = 299792458;   // m/s
      return mom / fourMom.E() * c; // m/s
   }

   static double dist(const XYZVector &x, const XYZVector &z) { return std::sqrt((x - z).Mag2()); }

   /// @brief to propagate the state vector using the process model
   /// @param x state vector
   /// @param v process noise vector
   /// @param vecZ The next measurement point used to stop the propagation.
   /// @return propagated (unaugmented) state vector
   static kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> &x, const kf::Vector<DIM_V> &v, const kf::Vector<DIM_Z> &z)
   {
      std::cout << "Staring to run funcF" << std::endl;
      // TODO: This needs to be filled with an RK4 solver for the physics model
      kf::Vector<DIM_X> y{x};
      XYZVector measurement(z[0], z[1], z[2]); // Measurement point in mm

      double charge = x[6];
      double eLoss = v[0];
      double mass = 938.272; // Mass in MeV/c^2

      double mat_density = 0; // Density of the material in g/cm^3
      AtTools::AtELossTable dedxModel(mat_density);
      dedxModel.LoadSrimTable(
         "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt"); // Load the
                                                                                                         // SRIM table
                                                                                                         // for energy
                                                                                                         // loss
      double scalingFactor = 1.0;
      int iterations = 0;
      double calc_eLoss = 0;

      while (std::abs(calc_eLoss - eLoss) > 1e-3) {
         std::cout << "Running iteration " << iterations << " with scaling factor: " << scalingFactor
                   << " and energy loss: " << calc_eLoss << std::endl;

         if (iterations > 100) {
            // If we are not converging, we should probably throw an error.
            throw std::runtime_error("Energy loss did not converge after 100 iterations.");
         }

         // Variables needed in a single run of the RK4 solver. This section needs to be repeated
         // until the energy loss converges to the correct value.
         double h = 1e-10; // Timestep in s (100 ns to start)
         double lastApproach = std::numeric_limits<double>::max();
         bool approaching = true;
         iterations++;

         XYZVector pos(x[0], x[1], x[2]);
         XYZVector mom(x[3], x[4], x[5]);
         double KE_initial = std::sqrt(mom.Mag2() + mass * mass) - mass; // Kinetic energy in MeV

         while (true) {
            std::cout << "Position: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << std::endl;
            std::cout << "Momentum: " << mom.X() << ", " << mom.Y() << ", " << mom.Z() << std::endl;

            // Using timestep, propagate state forward one step.
            double KE = std::sqrt(mom.Mag2() + mass * mass) - mass; // Kinetic energy in MeV
            auto dedx = scalingFactor * dedxModel.GetdEdx(KE);      // Get the stopping power in MeV/mm
            // std::cout << "KE: " << KE << " dedx: " << dedx << std::endl;

            auto spline = dedxModel.GetSpline();
            // std::cout << "Spline: " << spline.get_x_min() << " to " << spline.get_x_max() << std::endl;
            // std::cout << "dxde " << spline(KE) << " dxde " << dedx << std::endl;

            auto x_k1 = GetVel(mom, mass);
            auto p_k1 = Force(pos, mom, charge, mass, dedx);
            // std::cout << "vel: " << x_k1 << " speed " << x_k1.R() << std::endl;
            // std::cout << "Force: " << p_k1 << std::endl;

            auto x_k2 = GetVel(mom + p_k1 * h / 2, mass);
            auto p_k2 = Force(pos + x_k1 * h / 2, mom + p_k1 * h / 2, charge, mass, dedx);
            // std::cout << "vel: " << x_k2 << " speed " << x_k2.R() << std::endl;
            // std::cout << "Force: " << p_k2 << std::endl;

            auto x_k3 = GetVel(mom + p_k2 * h / 2, mass);
            auto p_k3 = Force(pos + x_k2 * h / 2, mom + p_k2 * h / 2, charge, mass, dedx);
            // std::cout << "vel: " << x_k3 << " speed " << x_k3.R() << std::endl;
            // std::cout << "Force: " << p_k3 << std::endl;

            auto x_k4 = GetVel(mom + p_k3 * h, mass);
            auto p_k4 = Force(pos + x_k3 * h, mom + p_k3 * h, charge, mass, dedx);
            // std::cout << "vel: " << x_k4 << " speed " << x_k4.R() << std::endl;
            // std::cout << "Force: " << p_k4 << std::endl;

            auto mom_SItoMeV = 1.60218e-13 / 299792458; // Factor to convert momentum to MeV/c from kg*m/s
            auto F_SI = (p_k1 + 2 * p_k2 + 2 * p_k3 + p_k4) / 6;

            // Convert momentum to SI, update, and convert back to MeV/c
            auto mom_SI = mom * mom_SItoMeV;
            mom_SI += F_SI * h;
            mom = mom_SI / mom_SItoMeV; // Convert back to MeV/c

            // Convert position to SI, update, and convert back to mm
            auto pos_SI = pos / 1e3; // Convert mm to m
            pos_SI += (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4) * h / 6;
            pos = pos_SI * 1e3; // Convert back to mm

            std::cout << "Average force: " << F_SI << " N" << std::endl;
            std::cout << "Momentum: " << mom * mom_SItoMeV << " kg m/s" << std::endl;
            std::cout << "Delta x: " << (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4) / 6 << " m/s" << std::endl;

            auto approach = dist(pos, measurement);
            std::cout << "pos: " << pos << " measurement" << measurement << std::endl;
            std::cout << "Approach: " << approach << " last approach: " << lastApproach << std::endl;
            if (approach < lastApproach) {
               // We are still approaching the measurement point
               approaching = true;
               lastApproach = approach;

               continue;
            }

            bool reachedMeasurementPoint = (approaching && approach > lastApproach);
            bool particleStopped = std::sqrt(mom.Mag2() + mass * mass) - mass < 0.01;
            if (reachedMeasurementPoint || particleStopped) {
               // Last iteration we were still approaching the measurement point. Now we are further away
               // then before. We have probably reached the measurement point if things are well behaved.
               // I can think of cases where this will not be true. A better solution might be to run
               // tracking the point of closest approach until the distance between the current state and
               // the measurement point is larger than the distance between the last state and the measurement point.

               y[0] = pos.X();
               y[1] = pos.Y();
               y[2] = pos.Z();
               y[3] = mom.X();
               y[4] = mom.Y();
               y[5] = mom.Z();

               // Update the scaling factor
               double KE_final = std::sqrt(mom.Mag2() + mass * mass) - mass;
               calc_eLoss = KE_initial - KE_final; // Energy loss in MeV
               scalingFactor *= calc_eLoss / eLoss;
               std::cout << "------- End of RK4 interation ---------" << std::endl;
               std::cout << "Particle stopped: " << particleStopped << std::endl;
               std::cout << "Reached measurement point: " << reachedMeasurementPoint << std::endl;
               std::cout << "Last approach: " << lastApproach << " Current approach: " << approach << std::endl;
               std::cout << "Desired energy loss: " << eLoss << " MeV" << std::endl;
               std::cout << "Calculated energy loss: " << calc_eLoss << " MeV" << std::endl;
               std::cout << "New scaling factor: " << scalingFactor << std::endl;
               std::cout << "Final Position: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << std::endl;
               std::cout << "Final Momentum: " << mom.X() << ", " << mom.Y() << ", " << mom.Z() << std::endl;
               return y;
            }
         } // End of loop over RK4 integration
      }    // End loop over energy loss convergence

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

   const double mass_p = 938.272;           // Mass of proton in MeV/c^2
   const double charge_p = 1.602176634e-19; // Charge of proton
};

// Definition of static member variables
TrackFitterUKFPhysicsTest::XYZVector TrackFitterUKFPhysicsTest::fBField;
TrackFitterUKFPhysicsTest::XYZVector TrackFitterUKFPhysicsTest::fEField;

TEST_F(TrackFitterUKFPhysicsTest, PhysicsPrediction)
{
   // TODO: This needs to be filled with a proper physics test.
   kf::Vector<DIM_X> x; // Initial state vector

   kf::Matrix<DIM_X, DIM_X> P; // Initial state vector covariance matrix

   // Note: process noise is defined in the UKF class, so we don't need to set it here.

   kf::Vector<DIM_Z> z; // Measurement vector to be used in the correction step
   z << std::sqrt(5), std::atan2(1.f, 2.f);

   kf::Matrix<DIM_N, DIM_N> R; // Covariance matrix for the measurement noise

   m_ukf.vecX() = x;
   m_ukf.matP() = P;

   // m_ukf.setCovarianceR(R);
   // m_ukf.predictUKF(funcF, z);
   ASSERT_EQ(true, true);
}

TEST_F(TrackFitterUKFPhysicsTest, TestForceNoFields)
{
   XYZVector pos(0, 0, 0);       // Position in mm
   XYZVector mom(100, 0, 0);     // Momentum in MeV/c
   fBField = XYZVector(0, 0, 0); // B-field in tesla
   fEField = XYZVector(0, 0, 0); // E-field

   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   double dedx = 1;          // Stopping power in MeV/mm

   auto force = Force(pos, mom, charge, mass, dedx);

   ASSERT_NEAR(force.X(), -1.602e-10, FLOAT_EPSILON);
   ASSERT_NEAR(force.Y(), 0, FLOAT_EPSILON);
   ASSERT_NEAR(force.Z(), 0, FLOAT_EPSILON);

   mom = XYZVector(100, 0, 100); // Reset momentum
   force = Force(pos, mom, charge, mass, dedx);
   ASSERT_NEAR(force.X(), -1.602e-10 / std::sqrt(2), FLOAT_EPSILON);
   ASSERT_NEAR(force.Y(), 0, FLOAT_EPSILON);
   ASSERT_NEAR(force.Z(), -1.602e-10 / std::sqrt(2), FLOAT_EPSILON);
}

TEST_F(TrackFitterUKFPhysicsTest, TestForceEField)
{
   XYZVector pos(0, 0, 0);       // Position in mm
   XYZVector mom(100, 0, 0);     // Momentum in MeV/c
   fBField = XYZVector(0, 0, 1); // B-field in tesla
   fEField = XYZVector(0, 0, 0); // E-field in V/m

   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   double dedx = 0;          // Stopping power in MeV/mm

   auto force = Force(pos, mom, charge, mass, dedx);

   ASSERT_NEAR(force.X(), 0, FLOAT_EPSILON);
   ASSERT_NEAR(force.Y(), 0, FLOAT_EPSILON);
   ASSERT_NEAR(force.Z(), 1.121e-14, FLOAT_EPSILON);
}

TEST_F(TrackFitterUKFPhysicsTest, TestForceBField)
{
   XYZVector pos(0, 0, 0);       // Position in mm
   XYZVector mom(100, 0, 0);     // Momentum in MeV/c
   fBField = XYZVector(0, 0, 1); // B-field in tesla
   fEField = XYZVector(0, 0, 0); // E-field

   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   double dedx = 0;          // Stopping power in MeV/mm

   auto force = Force(pos, mom, charge, mass, dedx);

   ASSERT_NEAR(force.X(), 0, FLOAT_EPSILON);
   ASSERT_NEAR(force.Y(), -5.09e-12, FLOAT_EPSILON);
   ASSERT_NEAR(force.Z(), 0, FLOAT_EPSILON);
}

TEST_F(TrackFitterUKFPhysicsTest, TestPropagatorNoField)
{

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   fBField = XYZVector(0, 0, 0);                  // B-field in tesla
   fEField = XYZVector(0, 0, 0);                  // E-field

   kf::Vector<DIM_X> x; // Initial state vector
   x[0] = 0;
   x[1] = 0;
   x[2] = 0;
   x[3] = p; // p_x
   x[4] = 0;
   x[5] = 0;

   ASSERT_NEAR(x[3], 43.331, 1e-1); // Make sure momentum is calculated correctly

   kf::Vector<DIM_V> v; // Process noise vector
   v[0] = 1;            // Energy loss in MeV
   v[1] = 0.0;          // No process noise in this example

   kf::Vector<DIM_Z> z; // Measurement vector
   z[0] = 1e3;
   z[1] = 0;
   z[2] = 0;

   auto final = funcF(x, v, z); // Propagate the state vector using the process model

   // Check the final position is close to the stopping point from LISE
   ASSERT_NEAR(final[3], 0, 0.1);  // Final momentum in x-direction should be close to 0
   ASSERT_NEAR(final[0], 210, 10); // Final position in x-direction should be close to 210 mm

   KE = 0.5;
   E = KE + mass_p;
   p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   x[3] = p;                               // Reset momentum
   final = funcF(x, v, z);                 // Propagate the state vector using the

   ASSERT_NEAR(final[3], 0, 0.1);  // Final momentum in x-direction should be close to 0
   ASSERT_NEAR(final[0], 68.6, 5); // Final position in x
}
