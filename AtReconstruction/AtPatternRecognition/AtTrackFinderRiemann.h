/********************************************************************
 * AtTrackFinderRiemann                                             *
 *                                                                  *
 * Smart-seed RANSAC with 3D inlier check.                           *
 *                                                                  *
 *   - For every RANSAC iteration: seed1 = random unassigned hit,    *
 *     seed2/seed3 = picked from seed1's k nearest neighbors in 3D   *
 *     (so the 3 seeds are physically local — not three hits drawn   *
 *     uniformly from the whole event, which would frequently span   *
 *     multiple tracks and produce a "compromise" circle).
 *   - Circle through 3 seeds (analytic) + a z = a + b*phi line      *
 *     anchored on the two seeds with the largest |Δphi|.            *
 *   - Inlier criterion is **3D**: |xy circle residual| < fInlierDist *
 *     AND |z residual| < fZInlierDist.                              *
 *   - Kasa LSQ refit on inliers (xy circle only); z-line is refit by *
 *     LSQ on the surviving inlier set; final inlier collection uses  *
 *     the refined circle and refined z-line.
 *                                                                  *
 * The historical name and class-id (3) are preserved so AtPRAtask
 * wiring stays backward compatible. Earlier revisions used a Riemann-
 * paraboloid lifted-space plane fit (hence the name) and later a
 * direct circle RANSAC with multiple post-filters; both produced too
 * much hit loss when filters chained. The current design moves the
 * 3D awareness up to seed selection and inlier check, eliminating
 * the need for downstream connectivity / Z-coherence post-filters.
 ********************************************************************/

#ifndef AtTRACKFINDERRIEMANN_H
#define AtTRACKFINDERRIEMANN_H

#include "AtPRA.h"

#include <Rtypes.h>

#include <memory>

class AtEvent;
class AtPatternEvent;

namespace AtPATTERN {

class AtTrackFinderRiemann : public AtPRA {
private:
   double fInlierDist{5.0};       ///< Max |xy circle residual| to accept hit (mm)
   double fZInlierDist{15.0};     ///< Max |z - (a + b*phi)| to accept hit along the helix (mm)
   int fMinHitsPerTrack{10};      ///< Reject candidates below this
   int fMaxTracks{6};             ///< Cap iterative extraction
   int fMaxIterations{400};       ///< RANSAC iterations per pass
   unsigned int fSeed{12345};     ///< RNG seed for reproducibility
   int fK_NN{10};                 ///< Number of nearest neighbors used as seed candidates
   double fMaxPhiGap{60.0};       ///< Max consecutive phi gap (deg) around fitted circle center before splitting an arc
   bool fUseSelectAndMerge{true}; ///< Apply AT-TPC primary/fragment heuristic (false for annular geometries like PUMA)

public:
   AtTrackFinderRiemann();
   ~AtTrackFinderRiemann() override = default;

   std::unique_ptr<AtPatternEvent> FindTracks(AtEvent &event) override;

   void SetInlierDist(double d) { fInlierDist = d; }
   void SetZInlierDist(double d) { fZInlierDist = d; }
   void SetMinHitsPerTrack(int n) { fMinHitsPerTrack = n; }
   void SetMaxTracks(int n) { fMaxTracks = n; }
   void SetMaxIterations(int n) { fMaxIterations = n; }
   void SetSeed(unsigned int s) { fSeed = s; }
   void SetKNN(int k) { fK_NN = k; }
   void SetMaxPhiGap(double deg) { fMaxPhiGap = deg; }
   void SetUseSelectAndMerge(bool use) { fUseSelectAndMerge = use; }

   ClassDefOverride(AtTrackFinderRiemann, 1);
};

} // namespace AtPATTERN

#endif
