#include "AtPropagator.h"

#include "AtELossTable.h"
#include "AtKinematics.h"

#include <Math/Vector3D.h>

#include <gtest/gtest.h>
#include <iostream>
#include <memory>
using ROOT::Math::Plane3D;
using ROOT::Math::XYZPoint;
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
   XYZPoint pos(0, 0, 0);    // Position in mm
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
   XYZPoint pos(0, 0, 0);    // Position in mm
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
   XYZPoint pos(0, 0, 0);    // Position in mm
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

TEST(AtPropagatorTest, PropagateToPoint_StoppingNoField)
{
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(
      "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt");
   AtPropagator propagator(charge, mass, std::move(elossModel));
   AtRK4Stepper stepper;
   AtMeasurementPoint measurementPoint({1e3, 0, 0});

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   XYZPoint startPos(0, 0, 0);                    // Start position in mm
   XYZVector startMom(p, 0, 0);                   // Start momentum in MeV/c

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   propagator.PropagateToMeasurementSurface(measurementPoint, stepper);

   auto finalPos = propagator.GetPosition();
   auto finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 210, 10); // Final position in x-direction should be close to 210 mm
   ASSERT_NEAR(finalMom.X(), 0, 0.1);

   KE = 0.75;
   E = KE + mass_p;
   p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   startMom.SetXYZ(p, 0, 0);               // Reset momentum
   propagator.SetState(startPos, startMom);

   propagator.PropagateToMeasurementSurface(measurementPoint, stepper); // Propagate to range
   finalPos = propagator.GetPosition();
   finalMom = propagator.GetMomentum();
   ASSERT_NEAR(finalPos.X(), 130, 10); // Final position in x-direction should be close to 130 mm
   ASSERT_NEAR(finalMom.X(), 0, 0.1);  // Final momentum in x-direction should be close to 0
}

TEST(AtPropagatorTest, PropagateToPoint_NoField)
{
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(
      "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt");
   AtPropagator propagator(charge, mass, std::move(elossModel));
   AtRK4Stepper stepper;
   AtMeasurementPoint measurementPoint({10, 0, 0});

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   XYZPoint startPos(0, 0, 0);                    // Start position in mm
   XYZVector startMom(p, 0, 0);                   // Start momentum in MeV/c

   double eLoss = 0.0285;                                     // Expected energy loss in MeV (LISE)
   double E_fin = KE - eLoss + mass_p;                        // Expected final energy after loss
   double p_fin = std::sqrt(E_fin * E_fin - mass_p * mass_p); // Expected final momentum in MeV/c

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   XYZPoint targetPoint(10, 0, 0); // Target point to propagate to 10 mm
   // propagator.PropagateToPoint(targetPoint, stepper);
   propagator.PropagateToMeasurementSurface(measurementPoint, stepper);

   auto finalPos = propagator.GetPosition();
   auto finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 10, 1); // Final position in x-direction should be close to 10 mm
   ASSERT_NEAR(finalMom.X(), p_fin, 0.01);
}

TEST(AtPropagatorTest, PropagateToPlane_NoField)
{
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(
      "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt");
   AtPropagator propagator(charge, mass, std::move(elossModel));
   AtRK4Stepper stepper;

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   XYZPoint startPos(0, 0, 0);                    // Start position in mm
   XYZVector startMom(p, 0, 0);                   // Start momentum in MeV/c

   double eLoss = 0.0285;                                     // Expected energy loss in MeV in 10 mm (LISE)
   double E_fin = KE - eLoss + mass_p;                        // Expected final energy after loss
   double p_fin = std::sqrt(E_fin * E_fin - mass_p * mass_p); // Expected final momentum in MeV/c

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   XYZPoint planePoint(10, 10, 10);        // Target point to propagate to 10 mm
   XYZVector planeNormal(1, 0, 0);         // Normal vector of the plane in x-direction
   Plane3D plane(planeNormal, planePoint); // Create the plane
   AtMeasurementPlane measurementPlane(plane);
   propagator.PropagateToMeasurementSurface(measurementPlane, stepper);

   auto finalPos = propagator.GetPosition();
   auto finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 10, 1); // Final position in x-direction should be close to 10 mm
   ASSERT_NEAR(finalMom.X(), p_fin, 0.1);
}

TEST(AtPropagatorTest, PropagateToPlane_StoppingNoField)
{
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(
      "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt");
   AtPropagator propagator(charge, mass, std::move(elossModel));
   AtRK4Stepper stepper;

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   XYZPoint startPos(0, 0, 0);                    // Start position in mm
   XYZVector startMom(p, 0, 0);                   // Start momentum in MeV/c

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   XYZPoint planePoint(220, 0, 0);         // Target point to propagate to 215 mm
   XYZVector planeNormal(1, 0, 0);         // Normal vector of the plane in x-direction
   Plane3D plane(planeNormal, planePoint); // Create the plane
   AtMeasurementPlane measurementPlane(plane);
   propagator.PropagateToMeasurementSurface(measurementPlane, stepper);

   auto finalPos = propagator.GetPosition();
   auto finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 220, 1); // Final position in x-direction should be close to 215 mm
   ASSERT_NEAR(finalMom.X(), 0, 0.1);
}

TEST(AtPropagatorTest, PropagateToPointAdaptive_NoField)
{
   double charge = charge_p; // Charge in Coulombs
   double mass = mass_p;     // Mass in MeV/c^2
   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(
      "/home/adam/fair_install/ATTPCROOTv2/AtReconstruction/AtFitter/OpenKF/kalman_filter/HinH.txt");
   AtPropagator propagator(charge, mass, std::move(elossModel));
   AtRK4AdaptiveStepper stepper;
   AtMeasurementPoint measurementPoint({10, 0, 0});

   double KE = 1; // Kinetic energy in MeV
   double E = KE + mass_p;
   double p = std::sqrt(E * E - mass_p * mass_p); // Momentum in MeV/c
   XYZPoint startPos(0, 0, 0);                    // Start position in mm
   XYZVector startMom(p, 0, 0);                   // Start momentum in MeV/c

   double eLoss = 0.0285;                                     // Expected energy loss in MeV in 10 mm (LISE)
   double E_fin = KE - eLoss + mass_p;                        // Expected final energy after loss
   double p_fin = std::sqrt(E_fin * E_fin - mass_p * mass_p); // Expected final momentum in MeV/c

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field
   propagator.SetH(1);              // Set initial step size to 1 s

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   propagator.PropagateToMeasurementSurface(measurementPoint, stepper);

   auto finalPos = propagator.GetPosition();
   auto finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 10, 10 * 1e-3); // Final position in x-direction should be close to 10 mm
   ASSERT_NEAR(finalMom.X(), p_fin, 0.1);

   propagator.SetState(startPos, startMom);
   propagator.SetEField({0, 0, 0}); // No electric field
   propagator.SetBField({0, 0, 0}); // No magnetic field
   propagator.SetH(1e-6);           // Set initial step size to 1e-6 m

   ASSERT_NEAR(propagator.GetMomentum().X(), 43.331, 1e-1);

   propagator.PropagateToMeasurementSurface(measurementPoint, stepper);

   finalPos = propagator.GetPosition();
   finalMom = propagator.GetMomentum();

   ASSERT_NEAR(finalPos.X(), 10, 10 * 1e-3); // Final position in x-direction should be close to 10 mm
   ASSERT_NEAR(finalMom.X(), p_fin, 0.1);
}