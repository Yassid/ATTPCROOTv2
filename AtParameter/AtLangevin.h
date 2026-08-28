#ifndef ATLANGEVIN_H
#define ATLANGEVIN_H

#include <cmath>

/**
 * @brief Langevin drift velocity of the ionisation electrons, in the PAD frame.
 *
 * There is exactly one place this is computed, and both directions of the calculation
 * must use it:
 *
 *   - AtClusterizeTask (simulation, forward)   displaces the electrons by +v_xy * t
 *   - AtPSA::CalculateXCorr/YCorr (reconstruction, reverse)   subtracts the same
 *
 * They were previously separate copies. That is dangerous rather than merely untidy: the
 * simulation's only claim to validating the correction is that truth and reconstruction
 * disagree by nothing, and if the two copies drift apart the residual silently measures
 * the difference between them instead of testing anything.
 *
 * Physics: eq. (4) of the AT-TPC commissioning paper (Bradt et al., NIM A 875 (2017) 65),
 *
 *     v_D = v/(1+wt^2) [ Ehat + wt (Ehat x Bhat) + wt^2 (Ehat.Bhat) Bhat ]
 *
 * The paper specialises this with B in the y-z plane, i.e. it assumes a tilt azimuth of
 * 90 deg. That assumption does not hold for the Dec 2014 data, where the beam -- and so B
 * -- sits at about -162 deg in the pad plane, so the azimuth is kept general here and
 * comes from AtDigiPar's ThetaRot. Passing azimDeg = 90 reproduces the published form.
 *
 * The result is finally rotated by ThetaPad, because the Langevin solution is derived in
 * the field frame while the coordinates it gets applied to are pad coordinates. Omitting
 * that rotation was the original bug: the shear came out with roughly the right magnitude
 * pointing ~110 deg away from where it belonged.
 */
namespace AtTools {

struct DriftVector {
   double x{0.}; ///< cm/us, pad frame
   double y{0.}; ///< cm/us, pad frame
   double z{0.}; ///< cm/us, along the drift axis
   double omegaTau{0.};
};

/**
 * @param vDrift    scalar drift velocity [cm/us]
 * @param bField    magnetic field [T]
 * @param eField    electric field [V/m]
 * @param tiltDeg   angle between E and B, i.e. the detector tilt [deg]
 * @param azimDeg   azimuth of B in the pad plane [deg]; 90 = the paper's y-z-plane form
 * @param thetaPadDeg  pad-plane rotation relative to the field frame [deg]
 */
inline DriftVector LangevinDrift(double vDrift, double bField, double eField, double tiltDeg,
                                 double azimDeg, double thetaPadDeg)
{
   DriftVector v;
   if (bField == 0. || eField <= 0.) { // no field: the drift is simply along the axis
      v.z = vDrift;
      return v;
   }
   const double deg = M_PI / 180.0;
   const double t = tiltDeg * deg, p = azimDeg * deg, tp = thetaPadDeg * deg;
   const double ot = (bField / eField) * vDrift * 1E4; // vDrift cm/us -> m/s
   const double front = vDrift / (1.0 + ot * ot);
   const double st = std::sin(t), ct = std::cos(t);
   const double sp = std::sin(p), cp = std::cos(p);

   // field frame
   const double fx = front * (ot * (-st * sp) + ot * ot * ct * st * cp);
   const double fy = front * (ot * (st * cp) + ot * ot * ct * st * sp);

   // rotate into the pad frame
   v.x = fx * std::cos(tp) - fy * std::sin(tp);
   v.y = fx * std::sin(tp) + fy * std::cos(tp);
   v.z = front * (1.0 + ot * ot * ct * ct);
   v.omegaTau = ot;
   return v;
}

} // namespace AtTools

#endif // ATLANGEVIN_H
