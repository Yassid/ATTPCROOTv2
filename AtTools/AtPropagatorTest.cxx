#include "AtPropagator.h"

#include "AtKinematics.h"

#include <Math/Vector3D.h>

#include <gtest/gtest.h>
#include <memory>
using ROOT::Math::XYZVector;

using namespace AtTools;

const double mass_p = 938.272;           // Mass of proton in MeV/c^2
const double charge_p = 1.602176634e-19; // Charge of proton

class DummyELossModel : public AtELossModel {
public:
   double eLoss = 1;
   DummyELossModel() : AtELossModel(0) {}

   double GetdEdx(double /*KE*/) const override { return eLoss; }
   double GetRange(double /*energyIni*/, double /*energyFin = 0*/) const override { return 1.0; }
   double GetEnergyLoss(double /*energyIni*/, double /*distance*/) const override { return 1.0; }
   double GetEnergy(double /*energyIni*/, double /*distance*/) const override { return 1.0; }
};

TEST(AtPropagatorTest, ForceNoField)
{
   XYZVector pos(0, 0, 0);   // Position in mm
   XYZVector mom(100, 0, 0); // Momentum in MeV/c

   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   double dedx = 1;          // Stopping power in MeV/mm

   // Create a dummy energy loss model
   auto elossModel = std::make_unique<DummyELossModel>();
   AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 0});
   propagator.SetBField({0, 0, 0});

   auto force = propagator.Force(pos, mom);

   ASSERT_NEAR(force.X(), -1.602e-10, 1e-12);
   ASSERT_NEAR(force.Y(), 0, 1e-12);
   ASSERT_NEAR(force.Z(), 0, 1e-12);

   mom = XYZVector(100, 0, 100); // Reset momentum
   force = propagator.Force(pos, mom);
   ASSERT_NEAR(force.X(), -1.602e-10 / std::sqrt(2), 1e-12);
   ASSERT_NEAR(force.Y(), 0, 1e-12);
   ASSERT_NEAR(force.Z(), -1.602e-10 / std::sqrt(2), 1e-12);
}

TEST(AtPropagatorTest, ForceEField)
{
   XYZVector pos(0, 0, 0);   // Position in mm
   XYZVector mom(100, 0, 0); // Momentum in MeV/c
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2

   // Create a dummy energy loss model
   auto elossModel = std::make_unique<DummyELossModel>();
   elossModel->eLoss = 0; // No energy loss for this test
   AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 70000});
   propagator.SetBField({0, 0, 0});

   auto force = propagator.Force(pos, mom);

   ASSERT_NEAR(force.X(), 0, 1e-12);
   ASSERT_NEAR(force.Y(), 0, 1e-12);
   ASSERT_NEAR(force.Z(), 1.121e-14, 1e-15);
}

TEST(AtPropagatorTest, ForceBField)
{
   XYZVector pos(0, 0, 0);   // Position in mm
   XYZVector mom(100, 0, 0); // Momentum in MeV/c
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2

   // Create a dummy energy loss model
   auto elossModel = std::make_unique<DummyELossModel>();
   elossModel->eLoss = 0; // No energy loss for this test
   AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 0});
   propagator.SetBField({0, 0, 1});

   auto force = propagator.Force(pos, mom);

   ASSERT_NEAR(force.X(), 0, 1e-12);
   ASSERT_NEAR(force.Y(), -5.09e-12, 1e-13);
   ASSERT_NEAR(force.Z(), 0, 1e-12);
}