#ifndef ATGENFITTER_H
#define ATGENFITTER_H
/**
 * @brief Clean genfit track fitter (a1975-era), no legacy machinery.
 *
 * A from-scratch replacement for AtFITTER::AtGenfit that keeps ONLY the working
 * genfit core (KalmanFitterRefTrack + RKTrackRep + AtSpacepointMeasurement) and
 * drops everything that caused trouble:
 *   - NO track merging / reclustering / single-vertex pool building
 *   - NO duplicated/inconsistent vertex-end logic, NO experiment enums
 *   - clusters are ordered by the DRIFT COORDINATE z (monotonic along any helix, so
 *     curling/spiraling tracks sequence correctly); the VERTEX is the highest-z_digi
 *     end (= lowest z_lab), the deterministic AT-TPC convention used by PRA/the UKF.
 *
 * Frame conventions mirror the validated UKF path (AtFitterUKF): lab z = ZPadPlane
 * - z_digi (uniform, no per-direction special case), no x reflection, and a signed
 * B field (default -2.85 T for experimental a1975 handedness). Derives from the
 * MODERN EventFit::AtFitter base so it plugs into the maintained AtFitterTask.
 *
 * AtGenfit is kept untouched as a fallback.
 */
#include "AtFitter.h"

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <Rtypes.h>

#include <memory>
#include <string>

class AtTrack;
class AtFittedTrack;
class AtFitMetadata;
class AtRawEvent;
class AtEvent;
class TClonesArray;

namespace genfit {
class AbsKalmanFitter;
class AbsMeasurement;
class AtSpacepointMeasurement;
template <class hit_T, class measurement_T>
class MeasurementProducer;
template <class measurement_T>
class MeasurementFactory;
} // namespace genfit
class AtHitCluster;

namespace AtTools {
class AtSpyralPID;
class AtParticleID;
} // namespace AtTools

namespace EventFit {

class AtGenfitter : public AtFitter {
public:
   AtGenfitter(Double_t bFieldTesla = -2.85, Int_t pdg = 1000010020, Double_t massAmu = 2.0135532, Int_t Z = 1,
               std::string eLossFile = "", Bool_t noMatEffects = kTRUE, Int_t minIter = 2, Int_t maxIter = 5);
   ~AtGenfitter();

   void Init() override;

   void SetParticle(Int_t pdg, Double_t massAmu, Int_t Z)
   {
      fPDG = pdg;
      fMassAmu = massAmu;
      fZ = Z;
   }
   void SetBField(Double_t bTesla) { fBField = bTesla; }
   void SetZPadPlane(Double_t z) { fZPadPlane = z; }
   void SetMinClusters(Int_t n) { fMinClusters = n; }
   void SetIterations(Int_t mn, Int_t mx)
   {
      fMinIter = mn;
      fMaxIter = mx;
   }
   void SetVertexAxisMaxDist(Double_t mm) { fVertexAxisMaxDist = mm; }

   /// Energy-loss model for genfit's MaterialEffects. Default (false) keeps the
   /// legacy AT-TPC behaviour: force the SRIM-table parameterization
   /// (useEnergyLossParam + setEnergyLossFile), which is accurate near the Bragg
   /// peak for the heavy ions the AT-TPC normally tracks (p, d, fragments) — this
   /// is what every validated heavy-ion pipeline (16C+p, (p,d), (d,p)) relies on.
   /// When true, use genfit's NATIVE Bethe-Bloch instead (skip the SRIM path):
   /// correct for minimum-ionizing / relativistic particles such as pions and
   /// kaons, whose kinetic energy lies outside the ion SRIM tables (the SRIM path
   /// throws "kinetic energy out of parameterization boundaries" on them). Use
   /// this for the PUMA channels. Opt-in so heavy-ion fits are byte-for-byte
   /// unchanged.
   void SetEnergyLossBetheBloch(Bool_t on) { fUseBetheBloch = on; }

   /// Per-track charge sign from the PRA curvature (AtTrack::GetChargeSign, the
   /// same cross-product prior the UKF uses). When on, each track is fit with the
   /// signed PDG matching its PRA charge: +|fPDG| if the PRA sign is >0, else
   /// -|fPDG| (e.g. fPDG=211 fits pi+ tracks and -211 fits pi- tracks). Needed for
   /// mixed-charge final states (PUMA branch 8 = pi+/pi-). Tracks with an unknown
   /// PRA sign (0) fall back to fPDG. Default off (single-species heavy-ion fits
   /// keep using fPDG unchanged). SetChargeSignFlip inverts the mapping if the
   /// genfit lab-frame handedness turns out opposite to PRA's reco-frame
   /// convention — determine empirically against truth, like the UKF did.
   void SetUseTrackChargeSign(Bool_t on) { fUseTrackChargeSign = on; }
   void SetChargeSignFlip(Bool_t on) { fChargeSignFlip = on; }

   /// Sign of the drift-z -> lab-z transform: z_lab = ZPadPlane + sign * z_digi.
   /// Default -1 is the experimental AT-TPC convention (z_lab = ZPadPlane - z_digi,
   /// z_digi a drift-time coordinate measured DOWN from the pad plane). PUMA sim
   /// digitizes hits already in a lab-like POSITIVE z frame (hits near the +75 mm
   /// vertex), so the -1 flip would push them to negative z, OUTSIDE the gas volume
   /// -> genfit material navigation fails and most fits die. Set +1 (with
   /// ZPadPlane=0) to use the drift z as-is. Opt-in; heavy-ion pipelines keep -1.
   void SetZDriftSign(Double_t s) { fZDriftSign = s; }

   /// Internal re-clustering (parallel to the UKF). PRA can leave only ~2 clusters
   /// per track on PUMA's fine annular pad plane, which starves genfit (needs
   /// >= fMinClusters well-separated points for a curvature fit). When on, ArcWalk
   /// re-clusters the raw hits into ~targetClusters clusters on a LOCAL COPY of the
   /// track, so the shared AtPatternEvent and any parallel UKF task are untouched.
   /// Same AtTrackTransformer::ClusterizeArcWalk the UKF uses. Default off (legacy
   /// heavy-ion fits keep using PRA's clusters directly).
   void SetReclusterArcWalk(Bool_t on, Int_t targetClusters = 8, Int_t minHitsPerCluster = 3, Int_t kNN = 10)
   {
      fReclusterArcWalk = on;
      fReclusterTarget = targetClusters;
      fReclusterMinHits = minHitsPerCluster;
      fReclusterKNN = kNN;
   }

   /// Resistive-pad ring centroiding: re-cluster the raw hits by annular RING
   /// (charge-weighted centroid per ring, one sub-pad point per ring) on a local
   /// copy, instead of ArcWalk. Exploits the resistive charge sharing (SetChargeDispersion
   /// in AtPulse) to beat pad quantization. Takes precedence over SetReclusterArcWalk.
   /// nPadsPerRing = azimuthal subdivision of the PUMA map (e.g. 256). Default off.
   void SetUseRingClustering(Bool_t on, Int_t nPadsPerRing = 256)
   {
      fRingClustering = on;
      fNPadsPerRing = nPadsPerRing;
   }

   /// Skip tracks whose PRA circle radius exceeds this (mm) BEFORE fitting — a
   /// near-straight track (radius -> inf) has unphysical momentum and makes genfit's
   /// RK stepper diverge/segfault (especially with precise ring-clustered centroids).
   /// The guard runs before the fit so it prevents the crash rather than catching it.
   /// e.g. 1500 mm => p > 1.8 GeV for |q|=1 at B=4 T, well above a 0.4 GeV pion.
   /// Default 0 (off).
   void SetMaxSeedRadius(Double_t mm) { fMaxSeedRadiusMM = mm; }

   /// Back-extrapolate the fitted vertex-end state to the point of closest approach
   /// to the beam axis (the z-axis line through the origin), using genfit's
   /// RKTrackRep::extrapolateToLine. When on, the AtFittedTrack reports:
   ///   - initialPositionXtr / vertex  = position AT the POCA (extrapolated)
   ///   - Kinematics (KE,theta,phi)    = momentum AT the POCA (post-extrap, at-vertex)
   ///   - initialPosition              = fitted state at the first cluster
   ///   - KinematicsXtr                = momentum at the first cluster (pre-extrap)
   /// matching the UKF convention so vertices are directly comparable. Without it,
   /// the vertex is just the fitted first-cluster state (~inner-ring radius from the
   /// axis). Default off. Only applies when material effects / a valid geometry are
   /// present (the extrapolation navigates the field); it silently keeps the
   /// first-cluster state if the extrapolation throws.
   void SetBackExtrapPOCA(Bool_t on) { fBackExtrapPOCA = on; }
   /// Keep genfit material effects ON during the POCA back-extrapolation, so the
   /// state is propagated through the real geometry (e.g. Cu trap + Al cryostat
   /// shells between the gas and the vertex) and the at-vertex momentum RECOVERS the
   /// energy the particle lost reaching the gas. Falls back to the geometric (no-
   /// material) extrapolation if the material propagation throws. Default off.
   void SetBackExtrapMaterial(Bool_t on) { fBackExtrapMaterial = on; }

   /// Constant subtracted from the back-extrapolated vertex z (mm), matching the
   /// UKF's SetVertexZBias. Compensates the hardwired +8.6 mm AtPSA shaping-delay
   /// offset in the digi z-reconstruction (peakingTime * vDrift). Only applied when
   /// SetBackExtrapPOCA is on. Default 0 (no correction) so other pipelines are
   /// unaffected. PUMA uses 8.6.
   void SetVertexZBias(Double_t mm) { fVertexZBias = mm; }

   /// Order clusters (and pick the vertex end) by distance from the beam axis (xy
   /// radius, ascending) instead of by drift z. The default z-ordering assumes the
   /// track advances monotonically in z (beam-along-z AT-TPC tracks), but PUMA's
   /// in-plane pions have pz ~ 0 -> z is nearly constant -> the z sort is degenerate
   /// and picks the OUTER end as the vertex ~half the time, seeding the fit backwards
   /// (it then diverges and the POCA lands far off-axis: the vertex dr-tail). Ordering
   /// by xy radius puts the beam-axis end first — correct for tracks that fan OUT from
   /// a near-axis vertex. Mirrors the UKF's SetUseClusterDirSeed xy sort. Default off
   /// (beam-along-z pipelines keep z-ordering).
   void SetVertexByXYRadius(Bool_t on) { fVertexByXYRadius = on; }

   /// Opt-in fix for BACKWARD-track seed reversal. The default seed always takes the
   /// lowest-z_lab end as the vertex and points the momentum toward +z_lab, which is
   /// correct for forward tracks but REVERSES backward ones (genfit reflects them into
   /// the forward hemisphere -> wrong theta/KE). When on, tracks PRA-tagged backward
   /// (GeoTheta > 90 deg) are instead seeded from the highest-z_lab end with the
   /// momentum pointing into the track (-z_lab). Forward tracks are unchanged. Default
   /// off, so the validated forward (p,d)/(p,p) pipelines are byte-for-byte preserved.
   void SetBackwardSeedFix(Bool_t on) { fBackwardSeedFix = on; }

   /// Acceptance window on the fitted polar angle (deg). Tracks outside are DROPPED
   /// as unphysical (near-beam / backward). Default 10-170 deg.
   void SetThetaWindow(Double_t minDeg, Double_t maxDeg)
   {
      fThetaMinDeg = minDeg;
      fThetaMaxDeg = maxDeg;
   }

   /// Per-cluster measurement resolution (mm) fed to the genfit fit. Sets the
   /// measurement covariance; absorbs multiple scattering when material effects are
   /// off. ~2-3 mm keeps real tracks converging while still constraining momentum.
   void SetMeasSigma(Double_t mm) { fMeasSigmaMM = mm; }

   /// Diffusion model for the per-cluster measurement covariance. The transverse (x,y)
   /// and drift-z position variance grow with drift distance L as
   ///   sigma^2 = sigma0^2 + D^2 * L_cm
   /// with sigma0 = SetMeasSigma and the diffusion coefficients D in mm/sqrt(cm). This
   /// captures the detector physics that long-drift clusters are blurrier (electron
   /// diffusion) so the fit weights near-pad-plane clusters more. The defaults (0, 0)
   /// reduce EXACTLY to the flat SetMeasSigma covariance, preserving the validated
   /// production; tune per gas/field (physics-dependent, A/B test before adopting).
   void SetDiffusion(Double_t transvMM_per_sqrtcm, Double_t longMM_per_sqrtcm)
   {
      fDiffTransMM = transvMM_per_sqrtcm;
      fDiffLongMM = longMM_per_sqrtcm;
   }
   /// Looseness factor applied to the drift-z measurement variance at zero drift
   /// (default 2.0: z from drift time is a bit less certain than the pad-plane x,y).
   void SetZLongFactor(Double_t f) { fZLongFactor = f; }
   /// If >0, scale each cluster's measurement variance by clamp(qRef/Q, 0.25, 4):
   /// high-charge (more-electron) clusters get a tighter centroid. Default 0 (off).
   void SetChargeCovRef(Double_t qRef) { fChargeRefForCov = qRef; }

   /// Post-fit quality cut. Tracks failing it are KEPT but flagged (GetGoodFit()==false)
   /// so curling/blob tracks remain inspectable. A fit passes if it converged, has
   /// ndf>0 with chi2/ndf < chi2NdfMax, and keMin < KE < keMax.
   void SetQualityCut(Double_t chi2NdfMax, Double_t keMin, Double_t keMax)
   {
      fChi2NdfMax = chi2NdfMax;
      fKEMin = keMin;
      fKEMax = keMax;
   }

   /// Only fit tracks whose Spyral PID (sqrtdEdx, brho) falls inside the gate loaded
   /// from this AtParticleID JSON (e.g. pid/deuteron_band.json). The fitter computes
   /// the PID itself (same AtSpyralPID the AtPIDTask uses), so it needs no PID branch.
   /// Out-of-gate tracks are skipped before the (expensive) genfit fit.
   void SetPIDGate(const std::string &jsonFile)
   {
      fPidGateFile = jsonFile;
      fUsePIDGate = kTRUE;
   }

   /// Continuity merging: before fitting, merge AtTrack fragments that PRA split
   /// across breaks but belong to ONE physical trajectory — clusters are concatenated
   /// and the fit's drift-z sort re-sequences them. Two fragments merge when:
   ///  - their PRA circles are the SAME circle: centres within centerDistMM (default
   ///    15 mm — fragments of one helix share a centre to ~10 mm; two distinct tracks,
   ///    even of similar radius, sit ~20+ mm apart, validated on a1975) and radii
   ///    within radiusFrac;
   ///  - their nearest endpoints are within gapMM in 3D AND that closest junction is
   ///    NOT vertex-end-to-vertex-end (which is two tracks sharing the production
   ///    vertex and diverging, not a continuation).
   /// Chains (A-B-C) merge transitively. Operates on a COPY of the pattern tracks, so
   /// the shared AtPatternEvent is untouched and a parallel UKF task on the same event
   /// is unaffected. Default off. (A global chord-direction guard was tried and dropped:
   /// it is anti-correlated on real data — strongly-curved split tracks have diverging
   /// chords while two short parallel tracks align — the shared-centre test is the clean
   /// discriminator.)
   void SetMergeContinuity(Bool_t on, Double_t centerDistMM = 15.0, Double_t radiusFrac = 0.3, Double_t gapMM = 30.0)
   {
      fMergeContinuity = on;
      fMergeCenterDist = centerDistMM;
      fMergeRadiusFrac = radiusFrac;
      fMergeGapMM = gapMM;
   }

   /// Event-level entry point. Overridden to optionally run continuity merging on a
   /// copy of the pattern tracks before the per-track fit loop.
   void FitEvent(AtTrackingEvent *trackingEvent, AtPatternEvent *patternEvent, AtFitMetadata *fitMetadata = nullptr,
                 AtRawEvent *rawEvent = nullptr, AtEvent *event = nullptr) override;

protected:
   AtFittedTrack *GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata = nullptr, AtRawEvent *rawEvent = nullptr,
                                 AtEvent *event = nullptr) override;

   /// Merge continuous AtTrack fragments (see SetMergeContinuity). Returns the merged
   /// track list; a no-op returning the input unchanged when merging is disabled.
   std::vector<AtTrack> MergeContinuousTracks(std::vector<AtTrack> tracks) const;

private:
   Double_t fBField;        // signed, Tesla (kGauss internally)
   Int_t fPDG;
   Double_t fMassAmu;
   Int_t fZ;
   std::string fELossFile;
   Bool_t fNoMatEffects;
   Int_t fMinIter, fMaxIter;
   Double_t fZPadPlane{1000.0};
   Int_t fMinClusters{4};
   Double_t fVertexAxisMaxDist{150.0}; // keep tracks whose closest cluster is within this of the beam axis
   Double_t fChi2NdfMax{5.0};          // post-fit quality cut thresholds (flag only, never drops)
   Double_t fKEMin{0.5};
   Double_t fKEMax{60.0};
   Double_t fMeasSigmaMM{4.0};         // per-cluster measurement resolution (mm)
   Double_t fDiffTransMM{0.0};         // transverse diffusion (mm/sqrt(cm)); var += fDiffTrans^2 * L_cm
   Double_t fDiffLongMM{0.0};          // drift-z diffusion (mm/sqrt(cm)); var += fDiffLong^2 * L_cm
   Double_t fZLongFactor{2.0};         // drift-z variance looseness factor at zero drift
   Double_t fChargeRefForCov{0.0};     // if >0, scale cov by clamp(qRef/Q): more charge -> tighter
   Double_t fThetaMinDeg{10.0};         // fitted-theta acceptance window (deg); outside -> dropped
   Double_t fThetaMaxDeg{170.0};
   Int_t fTPCDetID{0};
   Bool_t fInit{kFALSE};

   Bool_t fBackwardSeedFix{kFALSE};    // seed backward (GeoTheta>90) tracks from the far end (see setter)
   Bool_t fUseBetheBloch{kFALSE};      // use genfit native Bethe-Bloch (pions/kaons) instead of SRIM tables
   Bool_t fUseTrackChargeSign{kFALSE}; // pick signed PDG per track from PRA charge sign (mixed q final states)
   Bool_t fChargeSignFlip{kFALSE};     // invert the PRA-sign -> PDG-sign mapping (empirical frame calibration)
   Double_t fZDriftSign{-1.0};         // z_lab = ZPadPlane + fZDriftSign * z_digi (-1 legacy exp, +1 sim frame)
   Bool_t fReclusterArcWalk{kFALSE};   // ArcWalk re-cluster raw hits on a copy before fitting (PUMA)
   Bool_t fRingClustering{kFALSE};     // resistive ring centroiding (charge-weighted per ring)
   Int_t fNPadsPerRing{256};           // azimuthal subdivision for ring clustering
   Double_t fMaxSeedRadiusMM{0.0};     // skip near-straight tracks (radius > this) before fit; 0=off
   Int_t fReclusterTarget{8};          // target clusters for the ArcWalk re-clustering
   Int_t fReclusterMinHits{3};         // min hits per cluster
   Int_t fReclusterKNN{10};            // kNN adjacency for ArcWalk
   Bool_t fBackExtrapPOCA{kFALSE};     // back-extrapolate fitted state to beam-axis POCA (UKF-comparable vertex)
   Bool_t fBackExtrapMaterial{kFALSE}; // keep material effects on in the back-extrap (recover Cu/Al energy loss)
   Double_t fVertexZBias{0.0};         // constant subtracted from back-extrap vertex z (mm); AtPSA shaping-delay fix
   Bool_t fVertexByXYRadius{kFALSE};   // order clusters / pick vertex by xy radius (in-plane tracks) not drift z

   Bool_t fMergeContinuity{kFALSE};    // merge PRA-split fragments of one track before fitting
   Double_t fMergeCenterDist{15.0};    // max PRA-circle centre distance to merge (mm) — "same circle" test
   Double_t fMergeRadiusFrac{0.3};     // max fractional PRA-radius difference to merge
   Double_t fMergeGapMM{30.0};         // max 3D endpoint gap to merge (mm)

   Bool_t fUsePIDGate{kFALSE};
   std::string fPidGateFile;
   std::unique_ptr<AtTools::AtSpyralPID> fSpyralPID;
   std::unique_ptr<AtTools::AtParticleID> fPidGate;

   std::shared_ptr<genfit::AbsKalmanFitter> fKalmanFitter;
   TClonesArray *fHitClusterArray{nullptr};
   TClonesArray *fGenfitTrackArray{nullptr};
   genfit::MeasurementProducer<AtHitCluster, genfit::AtSpacepointMeasurement> *fMeasurementProducer{nullptr};
   genfit::MeasurementFactory<genfit::AbsMeasurement> *fMeasurementFactory{nullptr};

   ClassDefOverride(AtGenfitter, 1);
};

} // namespace EventFit
#endif
