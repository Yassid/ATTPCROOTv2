#include "AtPUMAGenerator.h"

#include <FairLogger.h>
#include <FairPrimaryGenerator.h>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>

#include <cmath>

ClassImp(AtPUMAGenerator);

namespace {
// Particle masses in GeV/c^2 (PDG 2024).
constexpr Double_t kMassKplus = 0.493677;
constexpr Double_t kMassPiMinus = 0.13957039;

// PDG codes.
constexpr Int_t kPdgKplus = 321;
constexpr Int_t kPdgPiMinus = -211;
} // namespace

AtPUMAGenerator::AtPUMAGenerator()
{
   SetName("AtPUMAGenerator");
}

Bool_t AtPUMAGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   // ---- Vertex ----------------------------------------------------------
   // Gaussian (x,y) clipped to a cylinder of radius fTrapRadius; uniform z.
   Double_t x0 = 0., y0 = 0.;
   for (int tries = 0; tries < 1000; ++tries) {
      x0 = gRandom->Gaus(fMeanX, fSigmaX);
      y0 = gRandom->Gaus(fMeanY, fSigmaY);
      if (std::sqrt(x0 * x0 + y0 * y0) <= fTrapRadius)
         break;
   }
   const Double_t z0 = fMeanZ + gRandom->Uniform(-fDeltaZ, fDeltaZ); // mm

   // ---- p-bar momentum (isotropic direction, Gaussian magnitude) -------
   const Double_t pmag = gRandom->Gaus(fMeanMomentum, fSigmaMomentum); // GeV/c
   Double_t ux = 0., uy = 0., uz = 0.;
   {
      const Double_t u = gRandom->Uniform(0., 1.);
      const Double_t v = gRandom->Uniform(0., 1.);
      const Double_t theta = std::acos(2. * u - 1.);
      const Double_t phi = 2. * TMath::Pi() * v;
      const Double_t st = std::sin(theta);
      ux = st * std::cos(phi);
      uy = st * std::sin(phi);
      uz = std::cos(theta);
   }
   const Double_t pbarPx = pmag * ux;
   const Double_t pbarPy = pmag * uy;
   const Double_t pbarPz = pmag * uz;

   // ---- W = target(He-3, at rest) + projectile(p-bar) -------------------
   TLorentzVector target(0., 0., 0., fTargetMass);
   const Double_t pbarE = std::sqrt(pmag * pmag + fBeamMass * fBeamMass);
   TLorentzVector projectile(pbarPx, pbarPy, pbarPz, pbarE);
   TLorentzVector W = target + projectile;

   // ---- Phase-space decay: K+ K+ pi- hexaquark --------------------------
   constexpr int kMult = 4;
   Double_t masses[kMult] = {kMassKplus, kMassKplus, kMassPiMinus, fHexMass};
   if (W.M() < masses[0] + masses[1] + masses[2] + masses[3]) {
      LOG(warn) << "AtPUMAGenerator: W = " << W.M() << " GeV below decay-product threshold "
                << masses[0] + masses[1] + masses[2] + masses[3] << " GeV. Skipping event.";
      return kFALSE;
   }
   TGenPhaseSpace event;
   event.SetDecay(W, kMult, masses);
   const Double_t weight = event.Generate();
   (void)weight; // unused; phase-space weight if needed for biased sampling

   const Int_t pdgs[kMult] = {kPdgKplus, kPdgKplus, kPdgPiMinus, fHexPdg};

   // ---- Push to primary generator (positions in cm per FairRoot convention)
   const Double_t vx_cm = x0 * 0.1;
   const Double_t vy_cm = y0 * 0.1;
   const Double_t vz_cm = z0 * 0.1;
   for (int i = 0; i < kMult; ++i) {
      const TLorentzVector *p = event.GetDecay(i);
      // FairPrimaryGenerator::AddTrack expects momentum in GeV/c.
      primGen->AddTrack(pdgs[i], p->Px(), p->Py(), p->Pz(), vx_cm, vy_cm, vz_cm);
   }

   return kTRUE;
}
