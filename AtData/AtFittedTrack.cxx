#include "AtFittedTrack.h"
#include "AtTools/AtKinematics.h"

#include <Rtypes.h>



ClassImp(AtFittedTrack);

using XYZVector = ROOT::Math::XYZVector;

const std::tuple<Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t> AtFittedTrack::GetEnergyAngles()
{
   using namespace AtTools::Kinematics;
   double fEnergy = GetKineticEnergy(fMomentum);
   double fEnergyXtr = GetKineticEnergy(fMomentumXtr);
   double fEnergyPRA = GetKineticEnergy(fMomentumGeo);

   return std::forward_as_tuple(fEnergy, fEnergyXtr, fMomentum.Theta(), fMomentum.Phi(), fEnergyPRA, fMomentumGeo.Theta(), fMomentumGeo.Phi());
}

AtFitter::FitStats AtFittedTrack::GetStats()
{
   return fFitStats;
}

const std::tuple<Int_t, Float_t, Float_t, Float_t, std::string> AtFittedTrack::GetTrackProperties()
{
   return std::forward_as_tuple(fTrackPropGeo.charge, fTrackPropGeo.brho, fTrackPropGeo.eloss, fTrackPropGeo.dedx, std::to_string(fTrackProp.pdg));
}
