#include "AtPiGunGenerator.h"

#include <FairLogger.h>
#include <FairPrimaryGenerator.h>

#include <TMath.h>
#include <TRandom.h>

#include <cmath>

ClassImp(AtPiGunGenerator);

namespace {
constexpr Double_t kMassPi_GeV = 0.13957039;
constexpr Int_t kPdgPiPlus = 211;
constexpr Int_t kPdgPiMinus = -211;
} // namespace

Bool_t AtPiGunGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   // ---- Vertex (mm). Gaussian (x,y) clipped to fTrapRadius cylinder,
   // uniform z in [fMeanZ ± fDeltaZ]. Falls back to (fMeanX, fMeanY) when
   // sigma=0 to support point-source vertices.
   Double_t x0 = fMeanX;
   Double_t y0 = fMeanY;
   if (fSigmaX > 0 || fSigmaY > 0) {
      for (int tries = 0; tries < 1000; ++tries) {
         x0 = gRandom->Gaus(fMeanX, fSigmaX);
         y0 = gRandom->Gaus(fMeanY, fSigmaY);
         if (std::sqrt(x0 * x0 + y0 * y0) <= fTrapRadius)
            break;
      }
   }
   const Double_t z0 = fMeanZ + (fDeltaZ > 0 ? gRandom->Uniform(-fDeltaZ, fDeltaZ) : 0.);

   // ---- Multiplicity ---------------------------------------------------
   const Int_t lo = std::max(1, fMinMult);
   const Int_t hi = std::max(lo, fMaxMult);
   const Int_t N = lo + (Int_t)gRandom->Uniform(0., (Double_t)(hi - lo + 1));

   // ---- Per-pion: random sign, isotropic direction, uniform |p| ------
   const Double_t vx_cm = x0 * 0.1;
   const Double_t vy_cm = y0 * 0.1;
   const Double_t vz_cm = z0 * 0.1;

   for (Int_t i = 0; i < N; ++i) {
      const Int_t pdg = (fForcedCharge > 0)   ? kPdgPiPlus
                        : (fForcedCharge < 0) ? kPdgPiMinus
                                              : (gRandom->Integer(2) == 0 ? kPdgPiPlus : kPdgPiMinus);

      // Isotropic direction: cos(theta) uniform in [-1, 1], phi uniform.
      const Double_t cosT = 2. * gRandom->Uniform(0., 1.) - 1.;
      const Double_t sinT = std::sqrt(std::max(0., 1. - cosT * cosT));
      const Double_t phi = 2. * TMath::Pi() * gRandom->Uniform(0., 1.);

      // Momentum magnitude uniform in [fMinP_GeV, fMaxP_GeV].
      const Double_t p_GeV = gRandom->Uniform(fMinP_GeV, fMaxP_GeV);
      const Double_t px = p_GeV * sinT * std::cos(phi);
      const Double_t py = p_GeV * sinT * std::sin(phi);
      const Double_t pz = p_GeV * cosT;

      primGen->AddTrack(pdg, px, py, pz, vx_cm, vy_cm, vz_cm);
   }

   LOG(debug) << "AtPiGunGenerator: emitted " << N << " pions at vertex (" << vx_cm << ","
              << vy_cm << "," << vz_cm << ") cm";
   (void)kMassPi_GeV;
   return kTRUE;
}
