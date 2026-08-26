#ifndef ATGAINMATCHTASK_H
#define ATGAINMATCHTASK_H

#include <FairTask.h>

#include <Rtypes.h>
#include <TClonesArray.h>
#include <TString.h>

#include <map>
#include <string>

/**
 * @brief Per-run dE/dx gain matching of an AtPIDEvent branch.
 *
 * The micromegas/GET gain drifts from run to run, so the same particle deposits a
 * different measured charge in run 19 than in run 133. Left alone, that smears the
 * sqrt(dE/dx)-vs-Brho plane across a run set and no single PID gate can be drawn on it.
 * This task rescales the PID dE/dx of each event by a factor tabulated per run:
 *
 *     dEdx      *= f(run)
 *     sqrtdEdx  *= sqrt(f(run))
 *
 * It is the ATTPCROOT counterpart of the custom Spyral `GainMatchPhase`, and it consumes
 * a table derived from the SAME measured scaling factors, so both pipelines normalise the
 * gain identically and a disagreement downstream cannot be blamed on this step.
 *
 * ORDER MATTERS: add it AFTER AtPIDTask. It modifies that task's AtPIDEvent branch in
 * place (via the public Clear/Add API -- no framework class is changed), so it must run
 * once the estimators have filled it. Adding it before AtPIDTask scales an empty event
 * and does nothing, silently.
 *
 * @note dE and arclength are deliberately NOT scaled, matching Spyral, which multiplies
 * only dEdx and sqrt_dEdx. The consequence is that dEdx != dE/arclength after gain
 * matching. If you recompute dE/dx from those two members downstream you get the
 * UNMATCHED value back, which is a real trap -- SetScaleDE(kTRUE) scales dE too and
 * restores the identity, at the cost of no longer matching Spyral exactly.
 */
class AtGainMatchTask : public FairTask {
private:
   TClonesArray *fPIDArray{nullptr}; ///< AtPIDEvent branch, modified in place

   TString fPIDBranchName{"AtPIDEvent"};
   TString fTablePath{""};

   Int_t fRun{-1};             ///< run number to look up; MUST be set
   Double_t fFactor{1.0};      ///< resolved factor for fRun
   Bool_t fScaleClassic{kTRUE};///< also scale the legacy AtPIDEstimator results
   Bool_t fScaleDE{kFALSE};    ///< also scale dE (breaks Spyral parity, restores dEdx=dE/arclen)
   Bool_t fAllowMissing{kFALSE};///< if the run is absent from the table, continue with f=1

   std::map<int, double> fTable;         ///< run -> factor
   std::map<int, std::string> fSource;   ///< run -> "measured" | "spline"

   Long64_t fNScaled{0};
   Long64_t fNEvents{0};

public:
   AtGainMatchTask() = default;
   AtGainMatchTask(const char *tablePath, Int_t run) : fTablePath(tablePath), fRun(run) {}
   ~AtGainMatchTask() override = default;

   /// CSV written by macro/Unpack_HDF5/C15d/tools/make_gainmatch_table.py
   /// ("run,factor,source", '#' comments allowed).
   void SetTable(const char *path) { fTablePath = path; }
   void SetRun(Int_t run) { fRun = run; }
   /// Parse the run number out of a "run_0017"-style name. Returns -1 if there is no
   /// digit group, which is a failure rather than a default -- gain matching the wrong
   /// run is worse than not gain matching at all.
   static Int_t RunNumberFromName(const TString &name);

   void SetPIDBranch(TString n) { fPIDBranchName = n; }
   void SetScaleClassic(Bool_t v) { fScaleClassic = v; }
   void SetScaleDE(Bool_t v) { fScaleDE = v; }
   void SetAllowMissingRun(Bool_t v) { fAllowMissing = v; }

   Double_t GetFactor() const { return fFactor; }

   InitStatus Init() override;
   void Exec(Option_t *opt) override;
   void Finish() override;

   ClassDefOverride(AtGainMatchTask, 1);
};

#endif // ATGAINMATCHTASK_H
