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
class AtELossModel;
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

   /// Opt-in fix for BACKWARD-track seed reversal. The default seed always takes the
   /// lowest-z_lab end as the vertex and points the momentum toward +z_lab, which is
   /// correct for forward tracks but REVERSES backward ones (genfit reflects them into
   /// the forward hemisphere -> wrong theta/KE). When on, tracks PRA-tagged backward
   /// (GeoTheta > 90 deg) are instead seeded from the highest-z_lab end with the
   /// momentum pointing into the track (-z_lab). Forward tracks are unchanged. Default
   /// off, so the validated forward (p,d)/(p,p) pipelines are byte-for-byte preserved.
   void SetBackwardSeedFix(Bool_t on) { fBackwardSeedFix = on; }

   /// Back-extrapolate the vertex-end state to the beam axis before reading position and
   /// momentum off it. genfit's getFittedState() with no argument is the FIRST MEASUREMENT
   /// POINT, not the reaction vertex: between them lies unmeasured gas in which the ejectile
   /// already lost energy, and that loss is never restored. AtFitterUKF has done this since
   /// the start (AtFitterUKF.cxx, "Back-extrapolation to beam axis") and on a1975 16C(d,t)15C
   /// it is worth a great deal -- with the correct D2 density it opens the reconstructed
   /// g.s.-to-3.103 spacing from 2.919 to 3.192 (truth 3.103) and cuts sigma 0.252 -> 0.155.
   ///
   /// ONLY MEANINGFUL WITH MATERIAL EFFECTS ON. With setNoEffects(true) the extrapolation
   /// transports the state geometrically and leaves |p| untouched, so it moves the vertex
   /// and changes nothing about the energy. It is also only as good as the material: the
   /// medium in the loaded geometry has to be the actual target gas, not a stand-in.
   /// Default off, so every existing pipeline is unchanged.
   void SetBackExtrapToAxis(Bool_t on) { fBackExtrapToAxis = on; }

   /// Add back, by hand, the energy the ejectile lost between the reaction vertex and the first
   /// measured cluster. Meant to be used WITH SetBackExtrapToAxis and WITHOUT genfit material
   /// effects: the extrapolation supplies the gap length geometrically, and CATIMA supplies the
   /// dE/dx over it, so the fit keeps the stability and the 85%-good-fit rate it has with
   /// material off (material on drops that to 60%) while still recovering the missing energy.
   /// Deliberately approximate -- one iteration on the mean energy over the gap, no straggling,
   /// no correction along the measured part of the track.
   /// @param density  g/cm3 of the gas (D2 at 300 torr, 293 K is 6.61e-5)
   /// @param matA     mass number of the target gas. 2 for deuterium; passing 1 here computes
   ///                 the loss in ordinary hydrogen, which for the same mass density is twice
   ///                 the stopping and is the mistake this whole exercise came out of.
   void SetManualELoss(Double_t density, Int_t matA);

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

   /// When a material-effects fit throws, retry the track with material effects OFF so it
   /// still gets a result instead of vanishing (default ON — the historical behaviour).
   /// Such tracks come from a different model and are marked in AtFitTrackMetadata via
   /// GetMatEffects()==false / GetMatEffectsFallback()==true; filter on that before quoting
   /// a resolution. Turn this OFF for a clean material-effects-only production, where a
   /// failed track is dropped rather than silently downgraded.
   void SetMatEffectsFallback(Bool_t on) { fMatEffectsFallback = on; }

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
   Bool_t fBackExtrapToAxis{kFALSE};   // extrapolate the vertex-end state to the beam axis (see setter)
   std::unique_ptr<AtTools::AtELossModel> fManualELoss; // optional hand-applied dE/dx over the vertex gap
   Bool_t fMatEffectsFallback{kTRUE};  // retry a failed matFX fit without material effects (flagged; see setter)

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
