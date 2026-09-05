#include "AtClusterize.h"

#include "AtDigiPar.h"
#include "AtMCPoint.h"
#include "AtSimulatedPoint.h"

#include <FairLogger.h>

#include <TClonesArray.h>
#include <TMath.h>     // for Sqrt, Cos, Sin, TwoPi
#include <TMathBase.h> // for Abs
#include <TObject.h>   // for TObject
#include <TRandom.h>
#include <TString.h> // for operator!=, TString

#include <algorithm> // for max
#include <utility>   // for move

thread_local AtClusterize::XYZPoint AtClusterize::fPrevPoint;
thread_local int AtClusterize::fTrackID = 0;

void AtClusterize::GetParameters(const AtDigiPar *fPar)
{
   fEIonize = fPar->GetEIonize() / 1000000; // [MeV]
   fFano = fPar->GetFano();
   fVelDrift = fPar->GetDriftVelocity();   // [cm/us]
   fCoefT = fPar->GetCoefDiffusionTrans(); // [cm^2/us]
   fCoefL = fPar->GetCoefDiffusionLong();  // [cm^2/us]

   fDetPadPlane = fPar->GetZPadPlane(); //[mm]
   fReverseDrift = fPar->GetReverseDrift();

   LOG(info) << "  Ionization energy of gas: " << fEIonize << " MeV";
   LOG(info) << "  Fano factor of gas: " << fFano;
   LOG(info) << "  Drift velocity: " << fVelDrift;
   LOG(info) << "  Longitudal coefficient of diffusion: " << fCoefL;
   LOG(info) << "  Transverse coefficient of diffusion: " << fCoefT;
   LOG(info) << "  Position of the pad plane (Z): " << fDetPadPlane;
   LOG(info) << "  Drift direction: " << (fReverseDrift ? "REVERSED (beam enters through the pad plane)"
                                                        : "normal (beam enters through the window)");
}

void AtClusterize::FillTClonesArray(TClonesArray &array, std::vector<SimPointPtr> &vec)
{
   for (auto &point : vec) {
      auto size = array.GetEntriesFast();
      new (array[size]) AtSimulatedPoint(std::move(*point)); // NO LINT
   }
}

std::vector<AtClusterize::SimPointPtr> AtClusterize::ProcessEvent(const TClonesArray &fMCPointArray)
{
   std::vector<SimPointPtr> ret;
   for (int i = 0; i < fMCPointArray.GetEntries(); ++i) {
      auto mcPoint = dynamic_cast<AtMCPoint *>(fMCPointArray.At(i));

      for (auto &point : processPoint(*mcPoint, i))
         ret.push_back(std::move(point));
   }
   return ret;
}

std::vector<AtClusterize::SimPointPtr> AtClusterize::processPoint(AtMCPoint &mcPoint, int pointID)
{
   if (mcPoint.GetVolName() != "drift_volume") {
      LOG(info) << "Skipping point " << pointID << ". Not in drift volume.";
      return {};
   }

   auto trackID = mcPoint.GetTrackID();
   XYZPoint currentPoint = getCurrentPointLocation(mcPoint); // [mm, mm, us]

   // If it is a new track entering the volume or no energy was deposited
   // record its location and the new track ID
   if (mcPoint.GetEnergyLoss() == 0 || fTrackID != trackID) {
      fPrevPoint = currentPoint;
      fTrackID = mcPoint.GetTrackID();
      return {};
   }

   auto genElectrons = getNumberOfElectronsGenerated(mcPoint);
   XYZVector step;
   if (genElectrons > 0) {
      step = (currentPoint - fPrevPoint) / genElectrons;
   }

   auto sigTrans = getTransverseDiffusion(currentPoint.z());  // mm
   auto sigLong = getLongitudinalDiffusion(currentPoint.z()); // us

   std::vector<SimPointPtr> ret;

   // Need to loop through electrons
   for (int i = 0; i < genElectrons; ++i) {
      auto loc = applyDiffusion(currentPoint + i * step, sigTrans, sigLong);
      XYZVector locVec(loc.X(), loc.Y(), loc.Z());
      ret.push_back(std::make_unique<AtSimulatedPoint>(pointID, i, locVec));
      LOG(debug2) << loc << " from " << currentPoint + i * step;
   }

   fPrevPoint = currentPoint;
   return ret;
}

double AtClusterize::getLongitudinalDiffusion(double driftTime)
{
   auto sigInCm = TMath::Sqrt(fCoefL * 2 * driftTime);
   auto sigInUs = sigInCm / fVelDrift;
   return sigInUs;
}

// Takes drift time in us returns sigma in mm
double AtClusterize::getTransverseDiffusion(double driftTime)
{
   return 10. * TMath::Sqrt(fCoefT * 2 * driftTime);
}

uint64_t AtClusterize::getNumberOfElectronsGenerated(const AtMCPoint &mcPoint)
{
   auto energyLoss = mcPoint.GetEnergyLoss() * 1000.;
   auto meanElec = energyLoss / fEIonize;
   auto sigElec = TMath::Sqrt(fFano * meanElec);
   // Clamp at zero: an unguarded Gaus(small_mean, sigma) can sample negative,
   // and an implicit cast to uint64_t wraps to ~1.8e19 — the for-loop in
   // processPoint then allocates AtSimulatedPoints until OOM.
   double n = gRandom->Gaus(meanElec, sigElec);
   return n > 0 ? static_cast<uint64_t>(n) : 0;
}

AtClusterize::XYZPoint AtClusterize::getCurrentPointLocation(const AtMCPoint &mcPoint)
{
   // THE ONE LINE THE REVERSED DETECTOR CHANGES.
   //
   // mcPoint.GetZ() is measured along the beam from where it enters the drift volume (0 at the
   // window end of the geometry, ZPadPlane at the far end). The drift length is the distance from
   // the ionisation to the PAD PLANE, so which end the pad plane is on is the whole difference:
   //
   //   normal   : pad plane at the far end  -> drift = ZPadPlane - z_beam
   //   reversed : pad plane at the near end -> drift = z_beam
   //
   // Nothing else here moves. The point is stored as (x, y, driftTime) and AtPulse works purely
   // in drift time, so the digitisation downstream of this line is identical in both modes; what
   // changes is how much diffusion a given track picks up and which time bucket it lands in.
   auto zInCm = fReverseDrift ? mcPoint.GetZ() : (fDetPadPlane / 10. - mcPoint.GetZ());
   auto driftTime = TMath::Abs(zInCm) / fVelDrift; // us

   return {mcPoint.GetX() * 10., mcPoint.GetY() * 10., driftTime};
}

AtClusterize::XYZPoint
AtClusterize::applyDiffusion(const AtClusterize::XYZPoint &loc, double_t sigTrans, double sigLong)
{
   auto r = gRandom->Gaus(0, sigTrans);
   auto phi = gRandom->Uniform(0, TMath::TwoPi());
   auto dz = gRandom->Gaus(0, sigLong);

   return loc + XYZVector(r * TMath::Cos(phi), r * TMath::Sin(phi), dz);
}
