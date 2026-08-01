#ifndef ATSPYRALPID_H
#define ATSPYRALPID_H

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>

#include <Rtypes.h>

#include <vector>

class AtTrack;

namespace AtTools {

/// Faithful C++ port of Spyral's EstimateResult (the physics observables produced
/// by spyral.core.estimator.estimate_physics).
struct AtSpyralResult {
   double brho{0.0};       // magnetic rigidity, T*m  = B*radius*0.001/|sin(polar)|
   double dEdx{0.0};       // charge over first arc inner-pad segment / spline arclength
   double sqrtdEdx{0.0};   // sqrt(|dEdx|)
   double dE{0.0};         // integrated charge over the first-arc inner segment
   double arclength{0.0};  // spline arc length of that segment (mm)
   double polar{0.0};      // polar angle (rad) from linregress(rho vs z)
   double azimuthal{0.0};  // azimuthal angle (rad)
   double radius{0.0};     // fitted first-arc circle radius (mm)
   int direction{-1};      // 0 forward, 1 backward, -1 none
   int nPoints{0};         // points used (after the z-tie collapse, so < nClusters)
   int nClusters{0};       // clusters on the source track -- the usual PID quality cut. Stored so a
                           // gate can be drawn from the PID output alone, without reopening the
                           // (much larger) pattern file just to count clusters.
   int trackID{-1};        // source AtTrack::GetTrackID(). NOTE this is AtTrackFinderTC's CLUSTER
                           // LABEL, not the track's array position -- a lone track can have ID 2.
                           // Match on it; never assume ID == index.
   ROOT::Math::XYZPoint vertex; // reaction vertex (mm)
   ROOT::Math::XYZPoint center; // spiral center (mm)
   bool valid{false};
   // debug intermediates (for ground-truth comparison vs real Spyral)
   double dbgSlope{0.0};
   int dbgMax1{-1};   // first-arc max index (pass 1, [1:]-relative)
   int dbgMax2{-1};   // first-arc max index (pass 2, over all points)
   int dbgTestIdx{0}; // linregress window size
};

/**
 * @brief Natural cubic smoothing spline, the C++ equivalent of scipy's
 *        make_smoothing_spline(x, y, lam).
 *
 * Minimizes  sum_i (y_i - f(x_i))^2 + lam * integral f''(t)^2 dt  over natural
 * cubic splines (Green & Silverman / Reinsch). Solves
 *   (T + lam Q^T Q) gamma = Q^T y,   g = y - lam Q gamma
 * with a dense linear solve (cluster sizes are small). Requires strictly
 * increasing x; Fit() returns false otherwise (mirrors Spyral's spline failure).
 */
class AtSmoothingSpline {
public:
   bool Fit(const std::vector<double> &x, const std::vector<double> &y, double lam);
   double Eval(double t) const;        // smooth value at arbitrary t
   const std::vector<double> &FittedValues() const { return fG; } // g_i at the knots
   bool IsValid() const { return fValid; }

private:
   std::vector<double> fX;  // knots (strictly increasing)
   std::vector<double> fG;  // fitted (smoothed) values at knots
   std::vector<double> fY2; // second derivatives at knots (natural: ends = 0)
   bool fValid{false};
};

/**
 * @brief Faithful port of Spyral's particle-ID estimator (estimate_physics_pass).
 *
 * Differs fundamentally from AtPIDEstimator: it isolates the FIRST ARC (up to the
 * point of maximum distance from the vertex), fits a dedicated circle to that arc,
 * derives the polar angle from a linear regression of rho-vs-z, refines the vertex
 * to the circle point nearest the beam axis, integrates charge over the inner-pad
 * portion of the first arc, and computes the arc length from a 1000-point spline
 * line integral. This naturally tames spiraling / multi-turn tracks.
 */
class AtSpyralPID {
public:
   AtSpyralPID() = default;

   AtSpyralResult Estimate(AtTrack &track) const;

   void SetBField(double b) { fBField = b; }
   void SetBeamRegionRadius(double r) { fBeamRegionRadius = r; }
   void SetSmoothingFactor(double s) { fSmoothing = s; }
   void SetDirectionThreshold(double t) { fDirThreshold = t; }
   void SetMinPoints(int n) { fMinPoints = n; }
   void SetSmallPadRadius(double r) { fSmallPadRadius = r; }
   /// Point source: hits = the AtHit point cloud (faithful to Spyral's cluster.data,
   /// ~100 pts/track); clusters = AtPRA's merged hit-clusters (~30 pts, coarser).
   void SetUseHits(bool h) { fUseHits = h; }
   void SetUseClusters(bool c) { fUseHits = !c; }
   /// Use the integrated pulse charge (AtHit::GetTraceIntegral) instead of the peak
   /// amplitude (GetCharge) for dE/dx. Matches Spyral's dE/dx scale (~x4). Default off.
   void SetUseTraceIntegral(bool u) { fUseTraceIntegral = u; }

private:
   // analytic least-squares circle fit (Spyral geometry.circle.least_squares_circle)
   static bool LeastSquaresCircle(const std::vector<double> &x, const std::vector<double> &y, double &xc, double &yc,
                                  double &radius);
   // Spyral clusterize.get_direction: winding sense of azimuth vs z
   int Direction(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z) const;

   double fBField{2.85};           // Tesla
   double fBeamRegionRadius{25.0}; // mm
   double fSmoothing{100.0};       // smoothing-spline lambda (estimation)
   double fDirThreshold{0.5};      // direction winding-fraction threshold
   int fMinPoints{30};             // min_total_trajectory_points
   double fSmallPadRadius{152.0};  // mm, inner-pad boundary for dE integration
   bool fUseHits{true};            // true = AtHit point cloud (faithful), false = hit clusters
   bool fUseTraceIntegral{false};  // true = integrated charge (Spyral scale), false = peak amplitude
};

} // namespace AtTools

#endif // ATSPYRALPID_H
