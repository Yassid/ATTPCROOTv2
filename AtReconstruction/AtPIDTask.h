#ifndef ATPIDTASK_H
#define ATPIDTASK_H

#include "AtHelix3D.h"
#include "AtPIDCombine.h"
#include "AtPIDEstimator.h"
#include "AtSpyralPID.h"

#include <FairTask.h>

#include <Rtypes.h>
#include <TClonesArray.h>
#include <TString.h>

#include <map>

/**
 * @brief Runs the two PID estimators on every pattern-recognized track.
 *
 * Reads the AtPatternEvent branch, and for each AtTrack computes PID observables
 * with BOTH the legacy AtTools::AtPIDEstimator (charge/arclength) and the faithful
 * Spyral-style AtTools::AtSpyralPID (first-arc + spline arclength). Results are
 * written to an AtPIDEvent branch (parallel per-track vectors) so the two methods
 * can be compared with full statistics. The legacy estimator is kept untouched.
 * It also runs AtTools::AtHelix3D (3D estimator, default ON, SetRunHelix3D(false) to skip), passing it
 * the Spyral direction as one of its three direction votes.
 *
 * Optional, all default OFF (2026-11-02; the a1975 D2 recipe, see docs/subsystems/pid.md):
 *  - SetPadGainMap: per-pad gain correction of the hit charges before every estimator runs;
 *  - SetCombine: one combined PID per track (AtPIDEvent::GetCombined), Spyral with a cross-check against the 3D
 *    circle Brho, the tracks failing it re-added from AtHelix3D (AtTools::CombinePID).
 */
class AtPIDTask : public FairTask {
private:
   TClonesArray *fPatternEventArray{nullptr}; ///< input (AtPatternEvent)
   TClonesArray fPIDArray;                     ///< output (AtPIDEvent)

   TString fInputBranchName{"AtPatternEvent"};
   TString fOutputBranchName{"AtPIDEvent"};
   Bool_t fIsPersistence{kTRUE};

   Double_t fBField{2.85};        ///< Tesla
   Double_t fSmallPadRadius{152.0};

   AtTools::AtPIDEstimator fClassicEst;
   AtTools::AtSpyralPID fSpyralEst;
   AtTools::AtHelix3D fHelix3DEst;
   Bool_t fRunHelix3D{kTRUE};

   std::map<int, double> fPadGain; ///< pad -> relative gain (empty = no correction)
   Bool_t fCombine{kFALSE};
   AtTools::AtPIDCombineOpt fCombineOpt;

public:
   AtPIDTask();
   virtual ~AtPIDTask() = default;

   void SetBField(Double_t b);
   void SetSmallPadRadius(Double_t r);
   /// AtSpyralPID::fMinPoints defaults to 30 and, until 2026-08-25, AtPIDTask had no way to
   /// change it -- so every persisted AtPIDEvent was frozen at 30 and could not be re-cut
   /// downstream. That default is the documented cause of ~99.94 % of PID rejections (knee at
   /// 15-20). Whatever value a production persists is the plane its gates must be drawn on.
   void SetMinPoints(Int_t n) { fSpyralEst.SetMinPoints(n); }
   void SetRunHelix3D(Bool_t v) { fRunHelix3D = v; }
   /// Per-run dE/dx gain correction (AtTools::AtPIDGainTable): sqrtdEdx of the Spyral and the 3D
   /// estimators multiplied by gSpyral and g3D. Default 1 (off). The legacy estimator is not touched.
   void SetDedxGain(Double_t gSpyral, Double_t g3D)
   {
      fSpyralEst.SetDedxGain(gSpyral);
      fHelix3DEst.SetDedxGain(g3D);
   }
   /// z_lab = zPadPlane - z_digi for the 3D estimator's lab-frame outputs (theta, vertex, backward)
   void SetZPadPlane(Double_t z) { fHelix3DEst.SetZTransform(-1, z); }
   /// direct access to tune the 3D estimator (blob removal, guards, expected rotation sense, ...)
   AtTools::AtHelix3D &GetHelix3D() { return fHelix3DEst; }
   /// 0 = exact-equality z-tie collapse (historical), which can make the arclength spline singular.
   void SetZTieTolerance(Double_t t) { fSpyralEst.SetZTieTolerance(t); }
   /// Use the integrated pulse charge instead of the peak amplitude for dE/dx, as Spyral does.
   /// ⚠ ONLY meaningful if the PSA actually measured the pulse width: AtPSAMultiFit must run with
   /// SetFloatPeakingTime(true), else the "integral" is a fixed multiple of the amplitude and this
   /// is a no-op up to an overall scale. The integral gives ~30-70 % more p/d/t/alpha separation
   /// but is NON-LINEAR in true dE/dx at high ionisation -- good for gating, not for dE/dx as a
   /// number. See macro/Unpack_HDF5/C15d/PID_ATTPCROOT_vs_SPYRAL.md.
   void SetUseTraceIntegral(Bool_t u) { fSpyralEst.SetUseTraceIntegral(u); }
   /// Divide every hit charge by its pad's relative gain before the estimators run (all three). The map is read with
   /// AtTools::AtPIDGainTable::LoadPadGains (`pad value...`, `column` = which value). Pads absent from the map are left
   /// alone. Returns false (and leaves the correction off) if the file cannot be read. Combine with the per-run gain
   /// (SetDedxGain): the pad map flattens the pad plane, the run factor follows the gain drift.
   Bool_t SetPadGainMap(TString file, Int_t column = 1);
   /// Fill AtPIDEvent::GetCombined: Spyral, unless its Brho fails the cross-check against the 3D circle, in which case
   /// the AtHelix3D result is used (see AtTools::AtPIDCombineOpt for the rule and the a1975 values). Needs the 3D
   /// estimator (SetRunHelix3D(true), the default); best with GetHelix3D().SetArcAngle(true) and SetHitRefit(true, true).
   void SetCombine(Bool_t on, const AtTools::AtPIDCombineOpt &opt = AtTools::AtPIDCombineOpt())
   {
      fCombine = on;
      fCombineOpt = opt;
   }
   void SetPersistence(Bool_t v) { fIsPersistence = v; }
   void SetInputBranch(TString n) { fInputBranchName = n; }
   void SetOutputBranch(TString n) { fOutputBranchName = n; }

   InitStatus Init() override;
   void Exec(Option_t *opt) override;

   ClassDefOverride(AtPIDTask, 2);
};

#endif // ATPIDTASK_H
