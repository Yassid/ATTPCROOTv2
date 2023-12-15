#ifndef ATFITTEDTRACK_H
#define ATFITTEDTRACK_H

#include <Math/Point3D.h>    // for PositionVector3D
#include <Math/Point3Dfwd.h> // for XYZPoint
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // for LorentzVector
#include <Rtypes.h>
#include <TObject.h>

#include <string> // for string
#include <tuple>  // for tuple

class TBuffer;
class TClass;
class TMemberInspector;
class AtTrack;

namespace AtFitter {
struct FitStats {

   Double_t pval;    //< P-Value of fit
   Double_t chi2;    //< Chi2 of fit in forward direction
   Double_t bchi2;   //< Chi2 of fit in backwards direction
   Double_t ndf;     //< Degrees of freedom of fit in forward direction
   Double_t bndf;    //< NDegrees of freedom of fit in backward direction
   Bool_t converged; //< Fit converged
};

struct TrackProp {
   Int_t pdg;            //< Particle ID in PDG code
   Double_t tracklength; //< Length of track from fit (mm)
};

struct TrackPropGeo {
   Int_t charge;   //< Charge state of particle
   Double_t brho;  //< Rigidity (nan if not defined)
   Double_t eloss; //< Average energy loss per hit
   Double_t dedx;  //< Estimation of average dedx from track
};

struct TrackPropExtr {
   Double_t poca;   //< Distance between extrapolated track and fInitialPosXtr
   Double_t length; //< Length of extrapolated track back between fInitialPos and fInitialPosXtr
};

} // namespace AtFitter

/**
 * This object represents the fit of a charged particle's track in the detector volume.
 * It contains both the physics information of the charge particle (p, E, and particle ID)
 * as well as the relative
 *
 */
class AtFittedTrack : public TObject {

   using XYZPoint = ROOT::Math::XYZPoint;
   using XYZEVector = ROOT::Math::PxPyPzEVector;

protected:
   Int_t fTrackID;      //< Track ID from pattern recognition
   Int_t fUniqueID{-1}; //< Unique ID. Set when added to AtTrackingEvent

   AtFitter::FitStats fFitStats;

   // Results of fitting method used
   XYZEVector fMomentum; //< 4-Momentum of track from fit (MeV)
   XYZPoint fInitialPos; //< Initial position of track in detector (mm)
   AtFitter::TrackProp fTrackProp;

   // Estimates of track parameters from geometry of track
   XYZEVector fMomentumGeo; //< 4-Momentum estimate from geometry of track (MeV)
   XYZPoint fInitialPosGeo; //< Initial position of track estimated from geometry (mm)
   AtFitter::TrackPropGeo fTrackPropGeo;

   XYZEVector fMomentumXtr; //< 4-Momentum estimate extrapolated to (0,0) (MeV)
   XYZPoint fInitialPosXtr; //< Initial position of track extrapolated to (0,0) (mm)
   AtFitter::TrackPropExtr fTrackPropExtr;
   /*
     Float_t fEnergy{0};
     Float_t fTheta{0};
     Float_t fPhi{0};
     Float_t fEnergyPRA{0};
     Float_t fThetaPRA{0};
     Float_t fPhiPRA{0};
     */

   /* Dependend on reaction
   Float_t fExcitationEnergy{0};
   Float_t fEnergyXtr{0};
   Float_t fExcitationEnergyXtr{0};
   */

public:
   AtFittedTrack();
   /**
    * @param[in] ID of the track (same as AtTrack object). Not necessarily unique.
    */
   AtFittedTrack(int trackID = -1){};
   AtFittedTrack(AtTrack &track){};

   ~AtFittedTrack() = default;

   inline void SetUniqueID(Int_t trackid) { fUniqueID = trackid; }
   /**
    * Sets the momentum, energy, and initial position of the fit track.
    */
   void SetFit(XYZEVector vec, XYZPoint point, AtFitter::TrackProp prop);

   /**
    * Sets the momentum, energy, and initial position of the track based on
    * geometry (typically from the pattern recognition with a circle).
    */
   void SetGeo(XYZEVector vec, XYZPoint point, AtFitter::TrackPropGeo prop);

   /**
    * Sets the momentum, energy, and initial position of the track extrapolating
    * the fit back to the point of closest approach.
    */
   void SetExtr(XYZEVector vec, XYZPoint point, AtFitter::TrackPropExtr prop);

   const Int_t GetTrackID() { return fTrackID; }
   const int GetID() { return fUniqueID; }
   AtFitter::FitStats GetStats();

   [[deprecated("Use getters for 4 momentum, 3 momentum, or energy")]] const std::tuple<
      Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t>
   GetEnergyAngles();

   [[deprecated("Use getters for track properties")]] const std::tuple<Int_t, Float_t, Float_t, Float_t, std::string>
   GetTrackProperties();

   ClassDef(AtFittedTrack, 2);
};

#endif