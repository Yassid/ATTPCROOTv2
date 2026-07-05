#ifndef ATPIDSTAMPTASK_H
#define ATPIDSTAMPTASK_H

#include "AtBetheBlochPID.h"

#include <FairTask.h>

#include <Rtypes.h>
#include <TString.h>

#include <string>
#include <vector>

class TClonesArray;

/**
 * @brief Stamps the authoritative species onto fitted tracks from a dE/dx vs rigidity PID.
 *
 * Runs AFTER the fitter(s). For each fitted track it looks up the source
 * pattern-recognized AtTrack (by TrackID), computes the magnetic rigidity and
 * truncated-mean dE/dx (AtTools::AtBetheBlochPID::Observables), classifies the species
 * with the Bethe-Bloch discriminant, and OVERWRITES the fitted track's ParticleInfo
 * (idPDG, charge sign preserved, mass) with the PID result. The charge sign is kept
 * from the fit/PRA (already curvature-determined); PID only fixes the mass identity —
 * exactly the information a same-rigidity kinematic fit cannot supply.
 *
 * The fitted-track branches are modified in place (no new branch); they are persisted
 * by the fitter task that registered them, so the stamped identity is written out.
 */
class AtPIDStampTask : public FairTask {
public:
   AtPIDStampTask();
   ~AtPIDStampTask() override = default;

   void SetBField(Double_t b) { fBField = b; }
   void SetPatternBranch(TString n) { fPatternBranch = n; }
   /// Fitted-track branch(es) to stamp (default UKF + genfit). Clears defaults on first call.
   void AddFittedBranch(const TString &n);

   /// Register a mass hypothesis (e.g. "pi",139.57039,211  and  "K",493.677,321).
   void AddSpecies(const std::string &name, double massMeV, int absPDG) { fPID.AddSpecies(name, massMeV, absPDG); }
   /// Convenience: register pi and K.
   void UsePionKaon();
   /// dE/dx-scale calibration (measured_dEdx = k * BetheBlochShape).
   void SetCalibration(double k) { fPID.SetCalibration(k); }
   void SetReferenceCalibration(double refRigidityMeV, double refDeDx, double refMassMeV)
   {
      fPID.SetReferenceCalibration(refRigidityMeV, refDeDx, refMassMeV);
   }
   /// Apply a pitch correction rigidity/sin(theta_fit) before classifying (default on).
   void SetPitchCorrection(Bool_t on) { fPitchCorrect = on; }

   InitStatus Init() override;
   void Exec(Option_t *opt) override;
   void Finish() override;

private:
   TClonesArray *fPatternArray{nullptr};
   std::vector<TClonesArray *> fFittedArrays;

   TString fPatternBranch{"AtPatternEvent"};
   std::vector<TString> fFittedBranches{"AtTrackingEventUKF", "AtTrackingEventGenfit"};
   bool fBranchesUserSet{false};

   Double_t fBField{4.0};
   Bool_t fPitchCorrect{kTRUE};
   AtTools::AtBetheBlochPID fPID;

   Long64_t fStamped{0};

   ClassDefOverride(AtPIDStampTask, 1);
};

#endif // ATPIDSTAMPTASK_H
