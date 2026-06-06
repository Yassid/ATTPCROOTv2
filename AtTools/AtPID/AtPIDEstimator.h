#ifndef ATPIDESTIMATOR_H
#define ATPIDESTIMATOR_H

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>

#include <Rtypes.h>

class AtTrack;

namespace AtTools {

/// PID observables for one track, mirroring spyral's EstimateResult columns.
struct AtPIDResult {
   double brho{0.0};      // magnetic rigidity, T*m  (B*R/sin(theta))
   double dEdx{0.0};      // integrated charge / arc length over the inner-pad segment
   double sqrtdEdx{0.0};  // sqrt(dEdx) — Spyral's default PID x-axis
   double elossMean{0.0}; // native pid_dev.C variable: mean cluster charge over 80% from the Bragg end
   double elossTrunc{0.0};// truncated mean (lowest 70%) of the same cluster charges
   double dE{0.0};        // integrated charge over the inner segment
   double arclength{0.0}; // arc length of the inner segment (mm)
   double polar{0.0};     // polar angle (rad, digi frame) from PRA
   double azimuthal{0.0}; // azimuthal angle (rad) from PRA
   int nClusters{0};      // number of clusters (or hits) used — for quality cuts
   ROOT::Math::XYZPoint vertex;
   bool valid{false};
};

/**
 * @brief Computes Spyral-style PID observables from a pattern-recognized AtTrack.
 *
 * C++ port of the relevant part of Spyral's estimator:
 *   brho = B * radius * 0.001 / |sin(polar)|            (radius from the PRA circle)
 *   dEdx = (charge integrated from the vertex outward,   (until r crosses the small-pad
 *           summed) / (arc length of that segment)        boundary, default r = 152 mm)
 *   sqrt_dEdx = sqrt(dEdx)
 *
 * The arc length is approximated by ordering the inner hits by cylindrical radius
 * and summing 3D segment lengths (Spyral uses a spline of the cluster arc; for the
 * near-vertex inner segment of recoil tracks the radial ordering is a close proxy).
 */
class AtPIDEstimator {
public:
   AtPIDEstimator(double bField_T = 2.85, double smallPadRadius_mm = 152.0)
      : fBField(bField_T), fSmallPadRadius(smallPadRadius_mm)
   {
   }

   AtPIDResult Estimate(AtTrack &track) const;

   void SetBField(double b) { fBField = b; }
   void SetSmallPadRadius(double r) { fSmallPadRadius = r; }
   double GetBField() const { return fBField; }
   double GetSmallPadRadius() const { return fSmallPadRadius; }

private:
   double fBField;        // Tesla
   double fSmallPadRadius; // mm, inner/outer pad boundary for dE/dx integration
};

} // namespace AtTools

#endif // ATPIDESTIMATOR_H
