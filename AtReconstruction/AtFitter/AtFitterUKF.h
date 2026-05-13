#ifndef ATFITTERUKF_H
#define ATFITTERUKF_H

#include "AtFitter.h"

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Rtypes.h>
#include <TMatrixD.h>

#include <memory>

class AtFitMetadata;
class AtFittedTrack;
class AtRawEvent;
class AtEvent;
class AtTrack;

namespace AtTools {
class AtELossModel;
}
namespace kf {
class TrackFitterUKF;
}

namespace EventFit {

/**
 * @brief FairTask-compatible UKF track fitter.
 *
 * Wraps kf::TrackFitterUKF and implements the EventFit::AtFitter interface so it
 * can be plugged into AtFitterTask.  For each AtTrack, the fitter:
 *   1. Seeds the momentum from the circle-fit Brho.
 *   2. Runs the UKF forward pass (predict+correct per cluster).
 *   3. Runs the RTS smoother.
 *   4. Returns an AtFittedTrack containing vertex kinematics and smoothed positions.
 *
 * Type convention (from CLAUDE.md): this class is never persisted to disk, so
 * plain C++ types (double/int/bool) are used — not ROOT typedefs.
 */
class AtFitterUKF : public AtFitter {
public:
   AtFitterUKF(double charge, double mass_MeV, std::unique_ptr<AtTools::AtELossModel> elossModel);
   ~AtFitterUKF(); // defined in .cxx so kf::TrackFitterUKF is complete there

   /// Set the magnetic field (Tesla).  May be called at any time; if the UKF is already initialized the update is propagated immediately.
   void SetBField(ROOT::Math::XYZVector bField);
   /// Set the electric field (V/m, default {0,0,0}).  May be called at any time; if the UKF is already initialized the update is propagated immediately.
   void SetEField(ROOT::Math::XYZVector eField);
   /// Set UKF sigma-point scaling parameters (alpha, beta, kappa).
   void SetUKFParameters(double alpha, double beta, double kappa);
   /// Set position measurement sigma (mm, default 1.0).
   void SetMeasurementSigma(double sigma_mm) { fMeasSigma_mm = sigma_mm; }
   /// Enable energy straggling in the propagator (default false).
   /// Only enable if the ELoss model supports GetRangeVariance() over the full energy range.
   void SetEnableEnergyStraggling(bool enable) { fEnableEnStraggling = enable; }
   /// Scale the forward-pass dE/dx in the propagator's process model
   /// (default 1.0 = no scaling). Diagnostic knob to test whether energy-loss
   /// model bias is responsible for residual KE bias; values <1 reduce dE/dx.
   void SetELossScaleFactor(double f) { fELossScaleFactor = f; }
   /// Set fractional momentum uncertainty for the initial covariance (default 0.1).
   void SetMomentumSigmaFrac(double frac) { fMomSigmaFrac = frac; }
   /// Minimum number of clusters required to attempt a fit (default 3).
   void SetMinClusters(int n) { fMinClusters = n; }
   /// Override the momentum seed (MeV/c). If > 0, this value is used instead of Brho.
   void SetMomentumSeed(double p_MeV) { fMomentumSeed = p_MeV; }

   /// Maximum back-extrapolation path length (mm) from the first smoothed
   /// cluster to the inferred vertex. The default (50 mm) is appropriate for
   /// beam+reaction runs where the vertex is close to the first cluster. For
   /// HELIOS-style runs with vertices deep inside the TPC, increase this
   /// (e.g. 200 mm) to let the linear extrapolation reach the true vertex.
   /// Set <= 0 to disable back-extrapolation entirely (vertex = first cluster).
   void SetBackExtrapMaxPath(double mm) { fBackExtrapMaxPath = mm; }
   /// Use the PRA circle (track->GetGeoCenter / GetGeoRadius) to compute the
   /// (x,y) POCA to the beam axis instead of the linear tangent step. Helix
   /// pitch carries the z; energy loss is integrated along the arc length.
   /// Required for setups with significant track curvature (e.g. PUMA at
   /// 4 T where K+ tracks make near-half loops). Default false preserves the
   /// 16C+p calibration where tracks are nearly straight and the propagator
   /// approach was reverted (commit aaabc31) due to Bragg-peak dE/dx
   /// instability in reverse — a problem that does NOT apply to PUMA's
   /// well-above-Bragg K+ and π−.
   void SetUseHelixBackExtrap(bool use) { fUseHelixBackExtrap = use; }
   /// Minimum spacing (mm) between adaptive clusters. Hits closer than this are
   /// merged. Default 3.0 mm matches AT-TPC pad pitch; lower for finer pad
   /// planes (e.g. PUMA's annular ~2 mm inner ring) where the default kills
   /// too many short, dense tracks via DROP[few-clusters].
   void SetMinClusterSpacing(double mm) { fMinClusterSpacing = mm; }
   /// Use the cluster-sequence derivative `clusters[1] - clusters[0]` as the
   /// direction seed for the UKF state, instead of PRA's GeoTheta/GeoPhi.
   /// Necessary for setups where PRA's (arclength,z) line fit ambiguates the
   /// forward direction (e.g. PUMA where `acos(sign·|Y|)` returns the
   /// supplementary angle for upward-going tracks). Default false preserves
   /// the legacy 16C+p path. Independent from `fZPadPlane`, which has its
   /// own cluster-derivative override active when set > 0.
   void SetUseClusterDirSeed(bool use) { fUseClusterDirSeed = use; }
   /// Switch the in-fitter re-clustering from `ClusterizeSmooth3D`
   /// (Z-seeded, gap-sensitive) to `ClusterizeArcWalk` (kNN spanning-tree
   /// arc ordering, geometry-based, gap-immune). Better for setups with
   /// centered vertices and arbitrary track directions (e.g. PUMA), where
   /// the highest-Z seed of Smooth3D lands at the wrong end of half the
   /// tracks. The target cluster count is taken from `fTargetClusters` (or
   /// the legacy 25 default if not set), `fArcWalkMinHits` is the floor on
   /// hits-per-cluster, and `fArcWalkKNN` controls graph adjacency.
   void SetUseArcWalk(bool use) { fUseArcWalk = use; }
   /// Force the back-extrapolated vertex to land at xy = (0, 0) instead of
   /// the POCA on the PRA circle. The PRA circle is fit through the hits
   /// only and typically does NOT pass through the actual beam-axis
   /// vertex (offset by a few mm), so the POCA-on-circle is a few mm off
   /// from the true vertex xy and the corresponding arc-length shortfall
   /// times cot(θ) shows up as a constant Δz bias (~+8 mm in PUMA).
   /// Enable for setups where the vertex is known to lie on the beam axis
   /// (e.g. PUMA's centered point source). 16C+p has vertex at high Z but
   /// also xy ≈ (0,0), so this is in principle safe there too — but it's
   /// kept as opt-in to avoid altering the calibrated 16C+p back-extrap
   /// behaviour. Default false.
   void SetForceVertexOnBeamAxis(bool force) { fForceVertexOnBeamAxis = force; }
   /// Rotate the smoothed φ by the back-extrapolation arc length divided by
   /// the PRA radius so AtFittedTrack::Kinematics.phi reports the angle AT
   /// THE VERTEX, consistent with KE at the vertex. By default the UKF
   /// stores the angle from the smoothed state at the FIRST CLUSTER, which
   /// inflates apparent σ_φ and bias_φ by the helix rotation between vertex
   /// and first cluster (event-by-event variable due to clustering wobble).
   /// θ is unchanged by Bz: motion in the cyclotron plane preserves the
   /// z-component, so only φ rotates. Sign uses fCharge: pi+/pi- handled
   /// automatically. Default false preserves all existing analyses.
   void SetUpdateAnglesOnBackExtrap(bool b) { fUpdateAnglesOnBackExtrap = b; }
   /// Subtract a constant offset (mm) from the back-extrapolated vertex z.
   /// Compensates a hardwired +peakingTime·vDrift bias in `AtPSA::CalculateZGeo`
   /// (the formula uses the pulse-peak time bucket but doesn't subtract the
   /// shaping-time delay that AtPulse added — these only cancel through
   /// `TBEntrance` in setups where the par file is calibrated for it).
   /// For PUMA: peakingTime=500 ns, vDrift=1.5 cm/μs ⇒ ~7.5 mm; an extra
   /// ~1 mm comes from pulse-shape asymmetry (measured constant +8.6 mm
   /// digi-vs-MC z offset across the full drift volume). Default 0.
   void SetVertexZBias(double mm) { fVertexZBias_mm = mm; }
   void SetArcWalkParams(int minHitsPerCluster, int kNN)
   {
      fArcWalkMinHits = minHitsPerCluster;
      fArcWalkKNN = kNN;
   }
   /// Override the default adaptive-clustering parameters (radius/distance, mm)
   /// passed to ClusterizeSmooth3D. Negative values keep the built-in
   /// energy-dependent defaults (25/15 below 3 MeV, 20/15 above), which were
   /// tuned for AT-TPC proton tracks. PUMA's short K+/π- tracks (~50-100 mm
   /// total) need tighter spacing (~5-7 mm) or all hits collapse into one
   /// cluster and hit DROP[few-clusters].
   void SetClusteringParams(double radius_mm, double distance_mm)
   {
      fClusterRadiusOverride = radius_mm;
      fClusterDistanceOverride = distance_mm;
   }
   /// Arc-length-aware adaptive clustering: target this many clusters per
   /// track regardless of length. Computes arc length from the PRA circle
   /// (R·|Δφ| spanned by hits around the GeoCenter), sets
   /// distance = clamp(arc/N, fAdaptiveDistMin, fAdaptiveDistMax) and
   /// radius = 1.5·distance. Set 0 to disable (use legacy KE-keyed defaults
   /// or the hard SetClusteringParams override). Tuned for PUMA where K+ and
   /// π- have wildly different track lengths (30-150 mm) but the legacy
   /// always-15 mm spacing collapses short tracks. Default 0.
   void SetTargetClusters(int n) { fTargetClusters = n; }
   /// Bounds on the per-track adaptive distance (mm). Defaults 4 mm / 20 mm
   /// keep clusters above sub-pad noise but below the AT-TPC default.
   void SetAdaptiveDistBounds(double minMm, double maxMm)
   {
      fAdaptiveDistMin = minMm;
      fAdaptiveDistMax = maxMm;
   }
   /// Number of fit iterations (default 1). With >1, the first pass uses wide covariance
   /// and subsequent passes use the previous result as seed with tighter covariance.
   void SetNIterations(int n) { fNIterations = n; }
   /// Use per-cluster covariance matrix from AtHitCluster::GetCovMatrix() instead of fixed sigma.
   void SetUsePerClusterCov(bool use) { fUsePerClusterCov = use; }
   /// Set the pad plane Z position (mm) for converting digi→lab coordinates.
   /// When set (>0), clusters are transformed to lab frame: Z_lab = ZPadPlane - Z_digi
   /// and sorted from vertex (high Z_lab) toward Bragg peak (low Z_lab).
   void SetZPadPlane(double z) { fZPadPlane = z; }
   /// Enable/disable adaptive re-clustering based on Brho momentum estimate.
   /// Disable when using pre-computed clusters (e.g. GNN arc-walk).
   void SetAdaptiveClustering(bool enable) { fAdaptiveClustering = enable; }

   /// AtFitter interface — no-op; UKF is created lazily on first GetFittedTrack() call.
   void Init() override {}

   /// Fit a single AtTrack and return a heap-allocated AtFittedTrack, or nullptr if
   /// the track has fewer than fMinClusters clusters or the UKF fails catastrophically.
   /// Public so multi-hypothesis wrappers (AtFitterUKFMulti) can call it directly.
   AtFittedTrack *GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata = nullptr, AtRawEvent *rawEvent = nullptr,
                                 AtEvent *event = nullptr) override;

protected:

private:
   /// Create and configure the kf::TrackFitterUKF.  Moves fELossModel to the propagator.
   /// Must only be called once.
   void InitUKF();

   ROOT::Math::XYZPoint GetInitialPosition(AtTrack *track) const;
   ROOT::Math::XYZVector GetInitialMomentum(AtTrack *track) const;
   TMatrixD GetInitialCovariance(double p_mag_MeV) const;
   TMatrixD GetMeasCovariance() const;
   double GetBrho(AtTrack *track) const;

   // Particle / detector configuration (plain C++ types — not persisted to disk)
   double fCharge;                            // Charge in Coulombs
   double fMass_MeV;                          // Mass in MeV/c^2
   ROOT::Math::XYZVector fBField{0, 0, 2.85}; // Magnetic field, Tesla
   ROOT::Math::XYZVector fEField{0, 0, 0};    // Electric field, V/m

   std::unique_ptr<AtTools::AtELossModel> fELossModel; // Moved to propagator on first use

   // UKF tuning parameters
   double fAlpha{1e-3};
   double fBeta{2.0};
   double fKappa{0.0};
   double fMeasSigma_mm{1.0};
   double fMomSigmaFrac{0.1};
   int fMinClusters{3};
   bool fEnableEnStraggling{true};
   double fELossScaleFactor{1.0};
   bool fUsePerClusterCov{false};
   double fZPadPlane{-1};      // If >0, convert digi→lab: Z_lab = fZPadPlane - Z_digi
   int fMaxFitTime_ms{2000};   // Maximum time per track fit in milliseconds
   double fMomentumSeed{-1};       // If >0, override Brho seed with this value (MeV/c)
   double fBackExtrapMaxPath{50.}; // mm; cap on back-extrapolation path. <=0 disables.
   bool fUseHelixBackExtrap{false}; // Circle-POCA xy + helix-pitch z (PUMA-style)
   double fClusterRadiusOverride{-1.0};   // mm; <=0 → use built-in defaults
   double fClusterDistanceOverride{-1.0}; // mm; <=0 → use built-in defaults
   int fTargetClusters{0};                // 0 → off; >0 enables arc-length-aware mode
   double fAdaptiveDistMin{4.0};
   double fAdaptiveDistMax{20.0};
   double fMinClusterSpacing{3.0}; // Minimum distance between clusters (mm)
   bool fUseClusterDirSeed{false}; // Override PRA Geo direction with cluster derivative
   bool fUseArcWalk{false};        // Use ClusterizeArcWalk instead of ClusterizeSmooth3D
   bool fForceVertexOnBeamAxis{false}; // Extrapolate to xy=(0,0), not POCA on PRA circle
   bool fUpdateAnglesOnBackExtrap{false}; // Rotate φ_s by arc/R so Kinematics.phi is at the vertex
   double fVertexZBias_mm{0.0};        // Subtract from final initialPositionXtr.Z()
   int fArcWalkMinHits{3};
   int fArcWalkKNN{10};
   double fMinLabTheta{10.0};      // Minimum lab angle (degrees) — skip beam-like tracks
   int fNIterations{1};            // Number of fit iterations (1=single pass)
   bool fAdaptiveClustering{true}; // Re-cluster based on Brho momentum estimate

   // Lazily initialised on first GetFittedTrack() call
   std::unique_ptr<kf::TrackFitterUKF> fUKF;
};

} // namespace EventFit

#endif // ATFITTERUKF_H
