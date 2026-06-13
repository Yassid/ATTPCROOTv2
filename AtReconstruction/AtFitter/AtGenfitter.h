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

protected:
   AtFittedTrack *GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata = nullptr, AtRawEvent *rawEvent = nullptr,
                                 AtEvent *event = nullptr) override;

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
   Double_t fThetaMinDeg{10.0};         // fitted-theta acceptance window (deg); outside -> dropped
   Double_t fThetaMaxDeg{170.0};
   Int_t fTPCDetID{0};
   Bool_t fInit{kFALSE};

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
