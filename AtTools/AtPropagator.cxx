#include "AtPropagator.h"

#include "AtKinematics.h"

#include <FairLogger.h>
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

void AtPropagator::RK4Step()
{
   double h = fH; // Step size in seconds

   auto x_k1 = Kinematics::GetVel(fMom, fMass);
   auto p_k1 = Force(fPos, fMom);

   auto x_k2 = Kinematics::GetVel(fMom + p_k1 * h / 2, fMass);
   auto p_k2 = Force(fPos + x_k1 * h / 2, fMom + p_k1 * h / 2);

   auto x_k3 = Kinematics::GetVel(fMom + p_k2 * h / 2, fMass);
   auto p_k3 = Force(fPos + x_k2 * h / 2, fMom + p_k2 * h / 2);

   auto x_k4 = Kinematics::GetVel(fMom + p_k3 * h, fMass);
   auto p_k4 = Force(fPos + x_k3 * h, fMom + p_k3 * h);

   auto F_SI = (p_k1 + 2 * p_k2 + 2 * p_k3 + p_k4) / 6; // Force in SI units (N)

   auto mom_SI = fReltoSImom * fMom;
   mom_SI += F_SI * h;          // Update momentum in SI units (kg m/s)
   fMom = mom_SI / fReltoSImom; // Convert back to

   auto pos_SI = fPos * 1e-3;                             // Convert position to SI units (m)
   pos_SI += (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4) * h / 6; // Update position in SI units (m
   fPos = pos_SI * 1e3;                                   // Convert back to mm
}

void AtPropagator::PropagateTo() {}

void AtPropagator::PropagateToPoint(const XYZPoint &point, double eLoss)
{
   LOG(info) << "Propagating to point: " << point << " with eLoss: " << eLoss;

   int iterations = 0;
   double calc_eLoss = 0;

   while (std::abs(calc_eLoss - eLoss) > 1e-4 || eLoss == 0) {
      LOG(debug) << "Running iteration " << iterations << " with scaling factor: " << fScalingFactor
                 << " and energy loss: " << calc_eLoss;

      if (iterations > 100) {
         // If we are not converging, we should probably throw an error.
         throw std::runtime_error("Energy loss did not converge after 100 iterations.");
      }

      double lastApproach = std::numeric_limits<double>::max();
      bool approaching = true;
      iterations++;
      auto KE_initial = Kinematics::KE(fMom, fMass);

      while (true) {
         auto lastPos = fPos;
         auto lastMom = fMom;
         LOG(debug) << "Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
         LOG(debug) << "Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();

         RK4Step();

         auto approach = (fPos - point).R();
         if (approach < lastApproach) {
            // We are still approaching the measurement point
            approaching = true;
            lastApproach = approach;
            continue;
         }

         bool reachedMeasurementPoint = (approaching && approach > lastApproach);
         bool particleStopped = Kinematics::KE(fMom, fMass) < fStopTol;
         if (reachedMeasurementPoint || particleStopped) {
            // Last iteration we were still approaching the measurement point. Now we are further away
            // then before. We have probably reached the measurement point if things are well behaved.
            // I can think of cases where this will not be true. A better solution might be to run
            // tracking the point of closest approach until the distance between the current state and
            // the measurement point is larger than the distance between the last state and the measurement point.

            // Undo the last step since we were closer last time.
            fPos = lastPos;
            fMom = lastMom;

            double KE_final = Kinematics::KE(fMom, fMass);
            calc_eLoss = KE_initial - KE_final; // Energy loss in MeV
            fScalingFactor *= eLoss / calc_eLoss;
            LOG(info) << "------- End of RK4 interation " << iterations << " ---------";
            LOG(info) << "Particle stopped: " << particleStopped;
            LOG(info) << "Reached measurement point: " << reachedMeasurementPoint;
            LOG(info) << "Last approach: " << lastApproach << " Current approach: " << approach;
            LOG(info) << "Desired energy loss: " << eLoss << " MeV";
            LOG(info) << "Calculated energy loss: " << calc_eLoss << " MeV";
            LOG(info) << "Difference: " << calc_eLoss - eLoss << " MeV";
            LOG(info) << "New scaling factor: " << fScalingFactor;
            LOG(info) << "Final Position: " << fPos.X() << ", " << fPos.Y() << ", " << fPos.Z();
            LOG(info) << "Final Momentum: " << fMom.X() << ", " << fMom.Y() << ", " << fMom.Z();

            if (eLoss == 0) {
               fScalingFactor = 1; // Reset scaling factor after convergence
               return;             // If no energy loss is specified, we are done.
            }
            break; // Else rerun with adjusted scaling factor
         }
      } // End of loop over RK4 integration
   }    // End loop over energy loss convergence

   fScalingFactor = 1; // Reset scaling factor after convergence
}
} // namespace AtTools
