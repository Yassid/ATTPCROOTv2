#include "AtPropagator.h"

#include "AtKinematics.h"
namespace AtTools {

AtPropagator::XYZVector AtPropagator::Force(XYZVector pos, XYZVector mom) const
{
   auto v = Kinematics::GetVel(mom, fMass);

   auto F_lorentz = fQ * (fEField + v.Cross(fBField));
   // std::cout << "F_lorentz: " << F_lorentz << std::endl;
   auto dedx = fScalingFactor * fELossModel->GetdEdx(Kinematics::KE(mom, fMass)); // Stopping power in MeV/mm
   auto dedx_si = dedx * 1.60218e-10;                                             // de_dx in SI units (J/m)

   auto drag = -dedx_si * mom.Unit();
   // std::cout << "drag: " << drag << " mom " << mom << " dedx " << dedx_si << std::endl;

   return F_lorentz + drag;
}
} // namespace AtTools
