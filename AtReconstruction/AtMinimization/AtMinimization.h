#ifndef ATMINIMIZATION_H
#define ATMINIMIZATION_H

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>

#include <memory>
#include <vector>

class AtEvent;
class AtHit;
class AtTrack;

/**
 * @brief Monte Carlo minimizers that fit a single track by simulating it.
 *
 * The classes in this namespace were ported from the FairRootv18.00 branch
 * (AtReconstruction/AtMinimization) written by Y. Ayyad and W. Mittig. The physics is
 * unchanged; the interfaces were adapted to the current data model (AtHit holding an
 * XYZPoint, events holding unique_ptr hits, AtMap owning the pad plane geometry).
 */
namespace MCMinimization {

/**
 * @brief Common interface for the Monte Carlo track minimizers.
 *
 * A minimizer takes a seed (the result of the pattern recognition) and a list of hits,
 * and searches for the track parameters whose simulation best reproduces those hits.
 */
class AtMinimization {
public:
   using XYZPoint = ROOT::Math::XYZPoint;
   using HitVector = std::vector<std::unique_ptr<AtHit>>;

   /**
    * @brief Starting point of the minimization.
    *
    * This is what the pattern recognition provides. It replaces the anonymous
    * `Double_t parameter[8]` array of the original implementation; the mapping is given
    * for each member so old macros can be translated.
    */
   struct TrackSeed {
      XYZPoint fVertex{};   //< Vertex position in the pad plane frame [mm] (parameter[0..2])
      int fVertexTB{0};     //< Time bucket of the vertex (parameter[3])
      double fPhi{0};       //< Azimuthal angle of the track [rad] (parameter[4])
      double fRadius{0};    //< Radius of curvature [mm], or the range [mm] if B == 0 (parameter[5])
      double fTheta{0};     //< Polar angle of the track [rad] (parameter[6])
      int fNumExpPoints{0}; //< Number of experimental points of the track (parameter[7])
   };

   /**
    * @brief Result of the minimization.
    *
    * Same content as the `FitPar` struct of the original implementation, with the `s`
    * prefix dropped and the units documented.
    */
   struct FitPar {
      double fTheta{0};       //< Polar angle at the minimum [rad] (sThetaMin)
      double fPhi{0};         //< Azimuthal angle at the minimum [rad] (sPhiMin)
      double fEnergy{0};      //< Energy per nucleon at the start of the track [MeV/u].
                              //< When B == 0 the track is seeded from its range and this is
                              //< the total kinetic energy [MeV]. (sEnerMin)
      XYZPoint fPos{};        //< Starting point of the fitted track [cm] (sPosMin)
      double fRadius{0};      //< Radius of curvature [mm], or the range [mm] if B == 0 (sBrhoMin)
      double fB{0};           //< Magnetic field at the minimum [G] (sBMin)
      double fChi2{0};        //< Objective function at the minimum (sChi2Min)
      double fChi2Q{0};       //< Charge/centroid term of the objective function (sChi2Q)
      double fChi2Range{0};   //< Range term of the objective function (sChi2Range)
      XYZPoint fVertexPos{};  //< Vertex from the backward extrapolation [cm] (sVertexPos)
      double fVertexEner{0};  //< Kinetic energy at the vertex [MeV] (sVertexEner)
      double fMinDistAppr{0}; //< Distance of minimum approach to the beam axis [cm] (sMinDistAppr)
      int fNumMCPoint{0};     //< Number of simulated points used in the position chi2 (sNumMCPoint)
      double fNormChi2{0};    //< Normalized chi2 (sNormChi2)
   };

protected:
   FitPar fFitPar{};

public:
   AtMinimization() = default;
   AtMinimization(const AtMinimization &) = default;
   AtMinimization(AtMinimization &&) = default;
   AtMinimization &operator=(const AtMinimization &) = default;
   AtMinimization &operator=(AtMinimization &&) = default;
   virtual ~AtMinimization() = default;

   /// Fit the track described by the passed hits. Returns false if the fit could not run.
   virtual bool Minimize(const TrackSeed &seed, const HitVector &hits) = 0;
   /// Convenience overload fitting every hit of an event.
   bool Minimize(const TrackSeed &seed, const AtEvent &event);
   /// Convenience overload fitting the hits of a single track.
   bool Minimize(const TrackSeed &seed, const AtTrack &track);

   const FitPar &GetFitPar() const { return fFitPar; }

   /// Clear the result of the previous fit.
   virtual void ResetParameters() { fFitPar = FitPar(); }
};

} // namespace MCMinimization

#endif // ATMINIMIZATION_H
