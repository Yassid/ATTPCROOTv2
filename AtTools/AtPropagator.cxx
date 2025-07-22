#include "AtPropagator.h"

#include "AtKinematics.h"

#include <FairLogger.h>

// Butcher tableau (c, a_ij) and the two b vectors for Dormand–Prince 5(4)
// c1 = 0
// Butcher tableau coefficients for Dormand–Prince 5(4) method

static constexpr double c[7] = {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};

static constexpr double a[7][6] = {
   {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
   {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
   {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
   {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0},
   {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0},
   {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0},
   {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}};

// b (5th-order)
static constexpr double b[7] = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};

// b* (4th-order, “star”)
static constexpr double bs[7] = {5179.0 / 57600.0, 0.0,       7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
                                 187.0 / 2100.0,   1.0 / 40.0};

using ROOT::Math::Plane3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;
namespace AtTools {

AtPropagator::XYZVector AtPropagator::Force(XYZPoint pos, XYZVector mom) const
{
   auto v = Kinematics::GetVel(mom, fMass);

   auto F_lorentz = fQ * (fEField + v.Cross(fBField));
   LOG(debug) << "F_lorentz: " << F_lorentz;
   auto dedx = fScalingFactor * fELossModel->GetdEdx(Kinematics::KE(mom, fMass)); // Stopping power in MeV/mm
   auto dedx_si = dedx * 1.60218e-10;                                             // de_dx in SI units (J/m)

   auto drag = -dedx_si * mom.Unit();
   LOG(debug) << "drag: " << drag << " mom " << mom << " dedx " << dedx_si;

   return F_lorentz + drag; // Force in N
}

AtPropagator::XYZVector AtPropagator::dpds(const XYZPoint &pos, const XYZVector &mom) const
{
   // Calculate the force acting on the particle at the given position and momentum
   auto speed = Kinematics::GetSpeed(mom.R(), fMass); // Speed in m/s
   return Force(pos, mom) / speed;
}
AtPropagator::XYZVector AtPropagator::d2xds2(const XYZPoint &pos, const XYZVector &mom) const
{
   auto phat = mom.Unit();         // Unit vector in the direction of momentum
   auto p = mom.R();               // Magnitude of the momentum
   auto dpds_vec = dpds(pos, mom); // Derivative of momentum w.r.t. arc length

   return 1 / p * (dpds_vec - phat * (phat.Dot(dpds_vec))); // Second derivative of position w.r.t. arc length
}

bool AtPropagator::ReachedPOCA(const XYZPoint &point)
{
   // Here we need to check if we are getting closer or further away from the POCA.
   // We may walk right past it so need to look for a change in the sign of the derivative or
   // something like that.
   auto lastDeriv = (fLastPos - point).Dot(fLastMom.Unit()); // proportional missing constants
   auto currDeriv = (fPos - point).Dot(fMom.Unit());
   LOG(debug) << "Last Derivative: " << lastDeriv << ", Current Derivative: " << currDeriv;
   return lastDeriv * currDeriv <= 0;
}

bool AtPropagator::IntersectedPlane(const Plane3D &plane)
{
   // Check if the particle has crossed the plane this step.
   auto prevSign = plane.Distance(fLastPos) > 0 ? 1 : -1;
   auto currSign = plane.Distance(fPos) > 0 ? 1 : -1;
   return (prevSign != currSign);
}

void AtPropagator::PropagateToPlane(const Plane3D &plane, AtStepper &stepper)
{
   LOG(info) << "Propagating to plane: " << plane;

   auto KE_initial = Kinematics::KE(fMom, fMass);
   stepper.fDeriv = [this](const XYZPoint &pos, const XYZVector &mom) { return this->Derivatives(pos, mom); };

   while (true) {
      LOG(debug) << "Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
      LOG(debug) << "Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();

      auto result = stepper.Step(fH, fPos, fMom);
      if (!result.success) {
         LOG(error) << "Integration step failed, aborting propagation.";
         return; // Abort propagation if step failed
      }
      CopyFromState(result); // Copy the new state from the stepper

      bool reachedMeasurementPoint = IntersectedPlane(plane);
      bool particleStopped = Kinematics::KE(fMom, fMass) < fStopTol;
      bool momentumReversed = (fLastMom.Dot(fMom) < 0);

      if (reachedMeasurementPoint && !particleStopped && !momentumReversed) {
         // We reached the measurement point, so we should figure out how far we are from the measurement point
         LOG(info) << "------ Reached measurement point ------";
         double finalH = (fLastPos - fPos).R(); // Distance traveled in the last step
         double approach = std::abs(plane.Distance(fLastPos));

         LOG(info) << "Distance to plane: " << approach << " mm";
         LOG(info) << "Final step size: " << finalH << " mm";

         finalH = approach * 1e-3; // Convert to meters for the RK4 step
         result = stepper.Step(finalH, fLastPos, fLastMom);
         if (!result.success) {
            LOG(error) << "Failed to propagate to measurement point, aborting.";
            return; // Abort propagation if step failed
         }
         auto origH = fH;       // Save original step size
         CopyFromState(result); // Update position and momentum to the new state
         fH = origH;            // Restore original step size
      }

      if (particleStopped || momentumReversed) {
         // In this case the particle stopped before hitting the plane
         // we should throw a warning to let the user know that there wasn't
         // enough energy to reach the plane.
         LOG(warning) << "------ Particle stopped before intersecting plane ------";

         // Calculate how far to travel before stopping
         double KE_last = Kinematics::KE(fLastMom, fMass);
         double deltaE = KE_last - fStopTol;
         deltaE = std::max(deltaE, 0.0); // Ensure we don't have negative energy loss

         LOG(info) << "Last KE: " << KE_last << " MeV";
         LOG(info) << "Energy to loose to stop: " << deltaE << " MeV";
         double h_Stop = deltaE / fELossModel->GetdEdx(KE_last); // Distance to stop in mm
         LOG(info) << "Estimated distance to stop: " << h_Stop << " mm";

         result = stepper.Step(h_Stop * 1e-3, fLastPos, fLastMom);
         if (!result.success) {
            LOG(error) << "Failed to propagate to stopping point, aborting.";
            return; // Abort propagation if step failed
         }
         auto origH = fH;       // Save original step size
         CopyFromState(result); // Update position and momentum to the new state
         fH = origH;            // Restore original step size
         LOG(info) << "Propagated to stopping point: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         LOG(info) << "Energy after stopping: " << Kinematics::KE(fMom, fMass) << " MeV";

         while (!IntersectedPlane(plane)) {
            fScalingFactor = 0; // Turn off energy loss.

            // If we still haven't intersected the plane, we need to adjust the step size
            double h = std::abs(plane.Distance(fPos)); // Reduce step size so we hit the plane
            if (h <= fDistTol)
               break;
            LOG(info) << "Propagating to plane after stopping with step size: " << h << " mm";
            result = stepper.Step(h * 1e-3, fPos, fMom);
            if (!result.success) {
               LOG(error) << "Failed to propagate to plane after stopping, aborting.";
               return; // Abort propagation if step failed
            }
            CopyFromState(result); // Update position and momentum to the new state
            LOG(info) << "New position after adjusting step size: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         }
         fLastMom = fMom;
         fMom = XYZVector(0, 0, 0); // Set momentum to zero since we stopped
         reachedMeasurementPoint = true;
      }

      if (reachedMeasurementPoint || particleStopped || momentumReversed) {
         double distanceToPlane = std::abs(plane.Distance(fPos));

         double KE_final = Kinematics::KE(fMom, fMass);
         auto calc_eLoss = KE_initial - KE_final; // Energy loss in MeV
         LOG(info) << "------- End of RK4 interation  ---------";
         LOG(info) << "Particle stopped: " << particleStopped;
         LOG(info) << "Reached measurement point: " << reachedMeasurementPoint;
         LOG(info) << "Distance to plane: " << distanceToPlane << " mm";
         LOG(info) << "Calculated energy loss: " << calc_eLoss << " MeV";
         LOG(info) << "Scaling factor: " << fScalingFactor;
         LOG(info) << "Final Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();

         // Project the position onto the plane. Cannot use ProjectOnPlane since it is templated in such
         // a way that it can't separate our XYZPoint and its internal XYZPoint.
         double d = plane.Distance(fPos); // Distance from the point to the plane
         fPos = XYZPoint(fPos.X() - plane.A() * d, fPos.Y() - plane.B() * d, fPos.Z() - plane.C() * d);
         LOG(info) << "Projected Position on plane: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         LOG(info) << "Final Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();
         return;
      }
   } // End of loop over RK4 integration
}

void AtPropagator::PropagateToMeasurementSurface(const AtMeasurementSurface &surface, AtStepper &stepper)
{
   LOG(info) << "Propagating to measurement surface";

   auto KE_initial = Kinematics::KE(fMom, fMass);
   stepper.fDeriv = [this](const XYZPoint &pos, const XYZVector &mom) { return this->Derivatives(pos, mom); };

   while (true) {
      LOG(debug) << "Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
      LOG(debug) << "Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();

      auto result = stepper.Step(fH, fPos, fMom);
      if (!result.success) {
         LOG(error) << "Integration step failed, aborting propagation.";
         return; // Abort propagation if step failed
      }
      CopyFromState(result); // Copy the new state from the stepper

      bool reachedMeasurementPoint = surface.PassedSurface(result);
      bool particleStopped = Kinematics::KE(fMom, fMass) < fStopTol;
      bool momentumReversed = (fLastMom.Dot(fMom) < 0);

      if (reachedMeasurementPoint && !particleStopped && !momentumReversed) {
         // We reached the measurement surface, so we should figure out how far we are from the measurement point
         LOG(info) << "------ Reached measurement surface ------";
         double finalH = (fLastPos - fPos).R(); // Distance traveled in the last step
         double approach = surface.Distance(fLastPos);

         LOG(info) << "Distance to plane: " << approach << " mm";
         LOG(info) << "Final step size: " << finalH << " mm";

         finalH = approach * 1e-3; // Convert to meters for the RK4 step
         result = stepper.Step(finalH, fLastPos, fLastMom);
         if (!result.success) {
            LOG(error) << "Failed to propagate to measurement point, aborting.";
            return; // Abort propagation if step failed
         }
         auto origH = fH;       // Save original step size
         CopyFromState(result); // Update position and momentum to the new state
         fH = origH;            // Restore original step size
      }

      if (particleStopped || momentumReversed) {
         // In this case the particle stopped before hitting the plane
         // we should throw a warning to let the user know that there wasn't
         // enough energy to reach the surface.
         LOG(warning) << "------ Particle stopped before reaching measurement surface ------";

         // Calculate how far to travel before stopping
         double KE_last = Kinematics::KE(fLastMom, fMass);
         double deltaE = KE_last - fStopTol;
         deltaE = std::max(deltaE, 0.0); // Ensure we don't have negative energy loss

         LOG(info) << "Last KE: " << KE_last << " MeV";
         LOG(info) << "Energy to loose to stop: " << deltaE << " MeV";
         double h_Stop = deltaE / fELossModel->GetdEdx(KE_last); // Distance to stop in mm
         LOG(info) << "Estimated distance to stop: " << h_Stop << " mm";

         result = stepper.Step(h_Stop * 1e-3, fLastPos, fLastMom);
         if (!result.success) {
            LOG(error) << "Failed to propagate to stopping point, aborting.";
            return; // Abort propagation if step failed
         }
         auto origH = fH;       // Save original step size
         CopyFromState(result); // Update position and momentum to the new state
         fH = origH;            // Restore original step size
         LOG(info) << "Propagated to stopping point: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         LOG(info) << "Energy after stopping: " << Kinematics::KE(fMom, fMass) << " MeV";

         while (surface.fClipToSurface && !surface.PassedSurface(result)) {
            fScalingFactor = 0; // Turn off energy loss.

            // If we still haven't intersected the surface, we need to adjust the step size
            double h = surface.Distance(fPos); // Reduce step size so we hit the surface
            if (h <= fDistTol)
               break;
            LOG(info) << "Propagating to surface after stopping with step size: " << h << " mm";
            result = stepper.Step(h * 1e-3, fPos, fMom);
            if (!result.success) {
               LOG(error) << "Failed to propagate to surface after stopping, aborting.";
               return; // Abort propagation if step failed
            }
            CopyFromState(result); // Update position and momentum to the new state
            LOG(info) << "New position after adjusting step size: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         }
         fLastMom = fMom;
         fMom = XYZVector(0, 0, 0); // Set momentum to zero since we stopped
         reachedMeasurementPoint = true;
      }

      if (reachedMeasurementPoint || particleStopped || momentumReversed) {
         double distanceToSurface = surface.Distance(fPos);

         double KE_final = Kinematics::KE(fMom, fMass);
         auto calc_eLoss = KE_initial - KE_final; // Energy loss in MeV
         LOG(info) << "------- End of RK4 interation  ---------";
         LOG(info) << "Particle stopped: " << particleStopped;
         LOG(info) << "Reached measurement point: " << reachedMeasurementPoint;
         LOG(info) << "Distance to plane: " << distanceToSurface << " mm";
         LOG(info) << "Calculated energy loss: " << calc_eLoss << " MeV";
         LOG(info) << "Scaling factor: " << fScalingFactor;
         LOG(info) << "Final Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();

         fPos = surface.ProjectToSurface(fPos);
         LOG(info) << "Projected Position on plane: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         LOG(info) << "Final Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();
         return;
      }
   } // End of loop over RK4 integration
}

void AtPropagator::PropagateToMeasurementSurface(const AtMeasurementSurface &surface, double eLoss, AtStepper &stepper)
{
   LOG(info) << "Propagating to surface with eLoss: " << eLoss;

   if (eLoss == 0) {
      LOG(warn) << "No energy loss specified, propagating without energy loss adjustment.";
      PropagateToMeasurementSurface(surface, stepper);
      return;
   }

   int iterations = 0;
   double calc_eLoss = 0;
   double KE_initial = Kinematics::KE(fMom, fMass);
   auto initialMom = fMom; // Save initial momentum for energy loss calculation
   auto initialPos = fPos; // Save initial position for energy loss calculation

   while (std::abs(calc_eLoss - eLoss) > 1e-4) {
      fMom = initialMom; // Reset position and momentum to initial values for the next iteration
      fPos = initialPos;

      LOG(debug) << "Running iteration " << iterations << " with scaling factor: " << fScalingFactor
                 << " and energy loss: " << calc_eLoss;

      if (iterations > 100) {
         // If we are not converging, we should probably throw an error.
         throw std::runtime_error("Energy loss did not converge after 100 iterations.");
      }

      iterations++;
      PropagateToMeasurementSurface(surface, stepper); // Propagate without energy loss adjustment

      double KE_final = Kinematics::KE(fMom, fMass);
      calc_eLoss = KE_initial - KE_final; // Energy loss in MeV
      fScalingFactor *= eLoss / calc_eLoss;
      LOG(info) << "Desired energy loss: " << eLoss << " MeV";
      LOG(info) << "Calculated energy loss: " << calc_eLoss << " MeV";
      LOG(info) << "Difference: " << calc_eLoss - eLoss << " MeV";
      LOG(info) << "New scaling factor: " << fScalingFactor;
      LOG(info) << "Condition: " << (std::abs(calc_eLoss - eLoss) > 1e-4);

   } // End loop over energy loss convergence

   LOG(info) << "Energy loss converged after " << iterations << " iterations.";

   fScalingFactor = 1; // Reset scaling factor after convergence
}

AtStepper::StepResult AtRK4Stepper::Step(double h, const XYZPoint &fPos, const XYZVector &fMom) const
{
   // Take h to be the step size in m.
   StepResult result;
   result.lastPos = fPos;
   result.lastMom = fMom;
   result.h = h;
   result.success = true;

   LOG(debug) << "Starting RK4 step with initial position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
   LOG(debug) << "Initial momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();
   LOG(debug) << "Step size (h): " << h << " m";

   auto [x_k1, p_k1] =
      fDeriv(fPos, fMom); // The derivative of the position is then just the unit vector of the momentum.

   auto x_2 = fPos + x_k1 * h / 2;               // Position at the midpoint
   auto p_2 = fMom + p_k1 * h / 2 / fReltoSImom; // Momentum at the midpoint

   auto [x_k2, p_k2] = fDeriv(x_2, p_2); // Derivative at the midpoint

   auto x_3 = fPos + x_k2 * h / 2;               // Position at the second midpoint
   auto p_3 = fMom + p_k2 * h / 2 / fReltoSImom; // Momentum at the second midpoint
   auto [x_k3, p_k3] = fDeriv(x_3, p_3);

   auto x_4 = fPos + x_k3 * h;               // Position at the end of the step
   auto p_4 = fMom + p_k3 * h / fReltoSImom; // Momentum at the end of the step
   auto [x_k4, p_k4] = fDeriv(x_4, p_4);

   auto dpds_SI = (p_k1 + 2 * p_k2 + 2 * p_k3 + p_k4) / 6; // "Force" in SI units (N)
   auto dxds_SI = (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4) / 6; // Position derivative in SI units (m)

   LOG(debug) << "dp/ds (SI units): " << dpds_SI.X() << ", " << dpds_SI.Y() << ", " << dpds_SI.Z();
   LOG(debug) << "dx/ds (SI units): " << dxds_SI.X() << ", " << dxds_SI.Y() << ", " << dxds_SI.Z();

   auto mom_SI = fReltoSImom * fMom;
   mom_SI += dpds_SI * h;             // Update momentum in SI units (kg m/s)
   result.mom = mom_SI / fReltoSImom; // Convert back to

   auto pos_SI = fPos * 1e-3; // Convert position to SI units (m)
   pos_SI += dxds_SI * h;     // Update position in SI units (m
   result.pos = pos_SI * 1e3; // Convert back to mm

   return result;
}

AtStepper::StepResult AtRK4AdaptiveStepper::Step(double h, const XYZPoint &fPos, const XYZVector &fMom) const
{
   // Take h to be the step size in m.
   StepResult result;
   result.lastPos = fPos;
   result.lastMom = fMom;
   result.h = h;
   result.success = true;

   // Take h to be the step size in m.
   // Use DP5(4) method for adaptive step size control.

   double atol_pos = 1e-2; // Absolute tolerance for position (mm)
   double atol_mom = 1e-2; // Absolute tolerance for momentum (MeV/c)
   double rtol = 1e-6;     // Relative tolerance for both position and momentum

   auto x0_mm = fPos;
   auto p0 = fMom;
   LOG(info) << "Starting RK4 step with initial position: " << x0_mm.X() << ", " << x0_mm.Y() << ", " << x0_mm.Z();
   LOG(info) << "Initial momentum: " << p0.X() << ", " << p0.Y() << ", " << p0.Z();

   while (true) {
      auto x_SI = fPos * 1e-3;        // Convert position to SI units (m)
      auto p_SI = fReltoSImom * fMom; // Convert momentum to SI units (kg m/s)
      XYZVector kx[7];                // kx[i] will hold the position derivatives (unitless)
      XYZVector kp[7];                // kp[i] will hold the momentum derivatives (SI units)

      // anonymous lambda to calculate and store the kx and kp values. Input is SI units.
      auto calc_k = [&](const XYZPoint &x, const XYZVector &p, int i) {
         auto [k_x, k_p] = fDeriv(x * 1e-3, p / fReltoSImom);
         kx[i] = k_x; // Store the position derivative (unitless)
         kp[i] = k_p; // Store the momentum derivative (SI units)
      };

      // anonymous lambda to calculate the position and momentum at the i-th stage
      auto calc_xp = [&](int i) {
         XYZVector dx(0, 0, 0);
         XYZVector dp(0, 0, 0);
         for (int j = 0; j < i; ++j) {
            dx = dx + kx[j] * a[i][j];
            dp = dp + kp[j] * a[i][j];
         }
         XYZPoint x = x_SI + dx * h;
         XYZVector p = p_SI + dp * h;
         return std::make_pair(x, p);
      };

      // Calculate kx and kp for each stage
      // build stage 0
      calc_k(x_SI, p0, 0);

      // build stage 1
      auto [x1, p1] = calc_xp(1);
      calc_k(x1, p1, 1); // k1

      // build stage 2
      auto [x2, p2] = calc_xp(2);
      calc_k(x2, p2, 2); // k2

      // build stage 3
      auto [x3, p3] = calc_xp(3);
      calc_k(x3, p3, 3); // k3

      // build stage 4
      auto [x4, p4] = calc_xp(4);
      calc_k(x4, p4, 4); // k4

      // build stage 5
      auto [x5, p5] = calc_xp(5);
      calc_k(x5, p5, 5); // k5

      // build stage 6
      auto [x6, p6] = calc_xp(6);
      calc_k(x6, p6, 6); // k6

      // Calculate the new position and momentum using the 5th-order method
      XYZVector dx(0, 0, 0);
      XYZVector dp(0, 0, 0);
      for (int i = 0; i < 7; ++i) {
         dx = dx + kx[i] * b[i];
         dp = dp + kp[i] * b[i];
      }
      XYZPoint x_new_5 = x_SI + dx * h;  // New position in SI units (m)
      XYZVector p_new_5 = p_SI + dp * h; // New momentum in SI units (kg m/s)

      // Calculate the new position and momentum using the 4th-order method
      dx = XYZVector(0, 0, 0);
      dp = XYZVector(0, 0, 0);
      for (int i = 0; i < 7; ++i) {
         dx = dx + kx[i] * bs[i];
         dp = dp + kp[i] * bs[i];
      }
      XYZPoint x_new_4 = x_SI + dx * h;  // New position in SI units (m)
      XYZVector p_new_4 = p_SI + dp * h; // New momentum in SI units (kg m/s)

      auto x_4_mm = x_new_4 * 1e3;          // Convert back to mm
      auto p_4_MeV = p_new_4 / fReltoSImom; // Convert back to MeV/c
      auto x_5_mm = x_new_5 * 1e3;          // Convert back to mm
      auto p_5_MeV = p_new_5 / fReltoSImom; // Convert back to MeV/c
      LOG(info) << "New position (5th order): " << x_5_mm.X() << ", " << x_5_mm.Y() << ", " << x_5_mm.Z();
      LOG(info) << "New momentum (5th order): " << p_5_MeV.X() << ", " << p_5_MeV.Y() << ", " << p_5_MeV.Z();
      LOG(info) << "New position (4th order): " << x_4_mm.X() << ", " << x_4_mm.Y() << ", " << x_4_mm.Z();
      LOG(info) << "New momentum (4th order): " << p_4_MeV.X() << ", " << p_4_MeV.Y() << ", " << p_4_MeV.Z();

      // Convert back to mm and MeV/c
      XYZVector x_err = (x_5_mm - x_4_mm);   // Error in position (mm)
      XYZVector p_err = (p_5_MeV - p_4_MeV); // Error in momentum (MeV/c)

      // Calculate the overall error
      double ex = x_err.X() / (atol_pos + rtol * std::abs(x_5_mm.X()));
      double ey = x_err.Y() / (atol_pos + rtol * std::abs(x_5_mm.Y()));
      double ez = x_err.Z() / (atol_pos + rtol * std::abs(x_5_mm.Z()));

      double ep_x = p_err.X() / (atol_mom + rtol * std::abs(p_5_MeV.X()));
      double ep_y = p_err.Y() / (atol_mom + rtol * std::abs(p_5_MeV.Y()));
      double ep_z = p_err.Z() / (atol_mom + rtol * std::abs(p_5_MeV.Z()));

      // Combine errors (norm)
      double err = std::sqrt(ex * ex + ey * ey + ez * ez + ep_x * ep_x + ep_y * ep_y + ep_z * ep_z);

      double factor = std::pow(err, -1.0 / 5.0); // Adjust step size based on error
      factor = std::clamp(factor, 0.25, 4.0);    // Clamp factor to reasonable limits
      double hNew = h * factor;
      // We now know the local error at this point. Now we need to decide to accept the point or not.
      if (err <= 1.0) {
         // Accept the step
         result.pos = x_5_mm;  // Update position in mm
         result.mom = p_5_MeV; // Update momentum in MeV/c
         LOG(info) << "Accepted step with error: " << err;
         LOG(info) << "Step size: " << h << " m";
         LOG(info) << "New step size: " << hNew << " m";
         LOG(info) << "New Position: " << result.pos.X() << ", " << result.pos.Y() << ", " << result.pos.Z();
         LOG(info) << "New Momentum: " << result.mom.X() << ", " << result.mom.Y() << ", " << result.mom.Z();

         // Adjust the step size for the next iteration
         result.h = hNew;
         result.success = true; // Step accepted
         return result;
      } else {
         // Reject the step and reduce the step size
         LOG(info) << "Rejected step with error: " << err;
         LOG(info) << "Step size: " << h << " m";
         LOG(info) << "Reducing step size to: " << hNew << " m";

         result.h = hNew; // Reduce step size
         h = hNew;        // Update the step size for the next iteration
         if (result.h < 1e-6) {
            LOG(error) << "Step size too small, aborting propagation.";
            result.success = false;
            return result; // Abort propagation if step size is too small
         }
      }
   }
}

bool AtMeasurementPoint::PassedSurface(AtStepper::StepResult &result) const
{
   // Check if the particle has passed the measurement point
   auto lastDeriv = (fPoint - result.lastPos).Dot(result.lastMom.Unit());
   auto currDeriv = (fPoint - result.pos).Dot(result.mom.Unit());
   LOG(debug) << "Last Derivative: " << lastDeriv << ", Current Derivative: " << currDeriv;
   return lastDeriv * currDeriv <= 0;
}

bool AtMeasurementPlane::PassedSurface(AtStepper::StepResult &result) const
{
   // Check if the particle has crossed the plane this step.
   auto prevSign = fPlane.Distance(result.lastPos) > 0 ? 1 : -1;
   auto currSign = fPlane.Distance(result.pos) > 0 ? 1 : -1;
   return (prevSign != currSign);
}

} // namespace AtTools
