#ifndef ATHELIX3D_H
#define ATHELIX3D_H

#include <Math/Point3D.h>
#include <Rtypes.h>

#include <vector>

class AtTrack;

namespace AtTools {

/// Result of AtHelix3D::Estimate, persisted in AtPIDEvent (one per pattern track, parallel to the
/// Spyral results). Lab frame: z_lab = zSign * z_digi + zOffset (see AtHelix3D::SetZTransform).
struct AtHelix3DResult {
   Bool_t valid{false};
   /// 0 ok; 1 fewer than minHits hits; 2 no backbone / end fits; 3 no PID segment on the small pads
   Int_t failCode{0};
   Int_t trackID{-1};  ///< source AtTrack::GetTrackID() -- match on it, never on the array position
   Int_t nClusters{0}; ///< hit clusters on the source track (the usual PID quality cut)

   // --- direction (vertex end) ---
   Int_t vertexEnd{-1};      ///< 0 = backbone end A, 1 = end B
   Bool_t ambiguous{true};   ///< the 2-of-3 vote tied or nobody voted -> the direction is unreliable
   Int_t nFor{0}, nVotes{0}; ///< votes agreeing with the chosen end / votes cast (max 3)
   Int_t voteSense{-1}, voteSpyral{-1}, voteBragg{-1}; ///< each voter's end, -1 = abstained
   Double_t braggRatio{0};   ///< charge(last third)/charge(first third) on the small pads, A -> B
   Double_t zLabVertexEnd{0}, zLabOtherEnd{0}; ///< lab z of the two track ends [mm]
   Bool_t backward{false};   ///< the vertex end has the HIGHER z_lab (the track travels toward -z_lab)

   // --- kinematics / PID from the vertex end ---
   Double_t theta{-1};          ///< polar angle of travel, lab frame [rad]
   Double_t brho{0};            ///< [T m]
   Double_t radius{0};          ///< vertex-end helix radius [mm]
   Double_t dca{1e9};           ///< closest approach of that helix to the beam axis [mm]
   ROOT::Math::XYZPoint vertex; ///< helix extrapolated to its closest approach to the axis (lab) [mm]
   Double_t dEdx{0};            ///< charge per mm of helix path, first arc on the small pads
   Double_t sqrtdEdx{0};
   Double_t pidPath{0}; ///< helix path used for dE/dx [mm]

   // --- track quality ---
   Int_t nHits{0}, nDropped{0};   ///< hits in, hits removed as stray blobs
   Double_t spacing{0};           ///< backbone length / backbone nodes [mm]; noise clouds are > ~6
   Double_t length{0};            ///< backbone length [mm]
   Bool_t guarded{false};         ///< vertex-end circle replaced by the whole-track circle
   Bool_t axisConstrained{false}; ///< vertex-end circle refitted through the beam axis

   // --- arc-length angle (SetArcAngle, default off; 2026-11-02) ---
   Double_t thetaArc{-1}; ///< polar angle from z_lab vs arc length round the vertex-end circle [rad]; -1 = not computed
   Double_t brhoArc{0};   ///< B * radius / sin(thetaArc) [T m]; the recommended 3D Brho (unbiased angle)

   // --- stopping / range quality (always computed when the backbone exists; 2026-11-02) ---
   Bool_t contained{false};  ///< far (stopping) end inside the gas: see AtHelix3D::SetContainment
   Double_t rStop{-1};       ///< radius of the far-end centroid [mm]
   Double_t zLabStop{-1};    ///< lab z of the far-end centroid [mm]
   Double_t peakFrac{-1};    ///< position of the charge-profile maximum / track length, from the vertex end
   Double_t tailRatio{-1};   ///< mean charge after the peak (+20 mm) / peak charge
   Bool_t midPeak{false};    ///< the charge peaks well before the end and the rest is weak: the particle stopped
                             ///< and the candidate carries other hits (see AtHelix3D::SetMidPeak)
};

/// Direction-free 3D track estimator: vertex end, Brho, theta and dE/dx without sorting hits by z.
///
/// WHY. AtPRA (sign of one RANSAC slope in arc-vs-z, hit-array order) and AtSpyralPID (azimuth steps
/// along the z-sorted cloud) both lean on z, which carries no ordering information near 90 deg.
/// This estimator orders hits in 3D (minimum-spanning-tree backbone), fits a helix at EACH end and
/// decides the vertex end by a 2-of-3 majority of independent calls. Algorithm: AtHelix3DCore.h.
///
/// VALIDATED on a1975 16C+d data, 10 runs (2026-10-30), against Spyral with kinematic referees:
///  - direction 96-97 % vs 93-95 % (elastic-locus referee); 100 % on 16C(d,p) simulation truth;
///  - deuteron ln-mass half-width 0.115 vs 0.148, outliers 8.9 vs 15.5 %, deuteron gate +9 points
///    efficiency at <= 3 % proton contamination; protons tie (0.114 vs 0.113).
///
/// CALIBRATION THAT IS SETUP-SPECIFIC: the expected rotation sense (+1 = counter-clockwise in x-y
/// travelling inward from the vertex end) was MEASURED as +1 for a1975 data AND simulation at
/// 2.85 T despite their opposite bFieldSign. Re-measure it for another field/charge/geometry.
class AtHelix3D {
public:
   /// @param spyralDirection AtSpyralResult::direction (0 forward = vertex at the LOWEST z_digi,
   ///        1 backward) or -1 if unavailable; it is one of the three direction votes.
   AtHelix3DResult Estimate(const AtTrack &track, int spyralDirection = -1) const;

   /// STRAY HITS (2026-11-02): indices into track.GetHitArray() of the hits the blob removal of this estimator drops (see
   /// SetBlobRemoval), except those farther than guardRadius [mm] from the beam axis: the wider spacing of the large pads
   /// opens gaps that split the backbone, so there the "blob" is often the real outer part of the track (on a1975 simulation
   /// deuterons at Brho 0.30-0.35 T m it was, and removing it ruined 3/4 of those tracks). Empty if the track has too few
   /// hits or no backbone could be built. Used by AtPRA::SetRemoveStrayHits.
   std::vector<int> FindStrayHits(const AtTrack &track, double guardRadius = 152.0) const;

   void SetBField(double b) { fBField = b; }                         ///< |B| [T]
   void SetZTransform(double sign, double offset) { fZSign = sign; fZOffset = offset; } ///< z_lab = sign*z_digi + offset
   void SetExpectedSense(int s) { fExpSense = s; }                    ///< see class comment
   void SetMinHits(int n) { fMinHits = n; }
   void SetSmallPadRadius(double r) { fSmallPadRadius = r; }          ///< dE/dx and Bragg use rho <= r
   void SetBinWidth(double w) { fBinW = w; }                          ///< centreline bin [mm]
   /// Blob removal: cut tree edges > |f| x the track's median edge (f < 0, floor 4 mm) or > f mm (f > 0);
   /// 0 disables it. Small = fewer than max(minBlob, frac*n) hits. Never drops more than maxDrop*n.
   void SetBlobRemoval(double f, int minBlob = 8, double frac = 0.15, double maxDrop = 0.5)
   {
      fMaxEdge = f; fMinBlob = minBlob; fBlobFrac = frac; fMaxDropFrac = maxDrop;
   }
   void SetBlobOnCircleTolerance(double mm) { fBlobOffCircle = mm; } ///< keep small pieces on the track circle
   void SetBlobLinearityKeep(double l) { fBlobLinMin = l; }          ///< keep line-like small pieces (0 = off)
   void SetRadiusGuard(double frac) { fGuardFrac = frac; }           ///< 0 disables
   /// Vertex-on-axis constraint: an end starting within rStart of the axis whose free circle misses
   /// the axis by more than dcaMax is refitted through the axis. dcaMax <= 0 disables it.
   void SetAxisConstraint(double dcaMax, double rStart = 45.0) { fDcaMax = dcaMax > 0 ? dcaMax : 1e9; fRStartAxis = rStart; }
   void SetBraggVote(double thr, double rContained = 150.0) { fBraggThr = thr; fRContained = rContained; }
   /// Per-run dE/dx gain correction: sqrtdEdx is multiplied by g (dEdx by g^2). Default 1 (off). The
   /// Bragg vote uses a charge RATIO and is unaffected. See AtPIDGainTable.
   /// Cap the end-helix (Brho) window at a fraction of the track and/or a length [mm] from the fitted end
   void SetEndWindowCap(double maxFrac, double maxLen = 1e9) { fWinMaxFrac = maxFrac; fWinMaxLen = maxLen; }
   /// Refit the end circles on the raw hits of their windows instead of the centroids (see AtHelix3DCore.h)
   void SetHitRefit(bool b, bool alsoZ = false) { fHitRefit = b; fHitRefitZ = alsoZ; }
   void SetDedxGain(double g) { fDedxGain = g; }
   double GetDedxGain() const { return fDedxGain; }
   /// ARC-LENGTH POLAR ANGLE (default off). theta from a straight-line fit of z_lab against the arc length round the
   /// vertex-end circle, over the first `frac` of the vertex-end window hits (at least 6), and brhoArc = B R/sin(theta).
   /// On a1975 simulation the angle is unbiased to 0.2-0.3 % in sin(theta), where the whole-window z vs phase angle
   /// biases slow tracks. Best with SetHitRefit(true, true). Needs at least 12 window hits.
   void SetArcAngle(bool on, double frac = 0.5) { fArcAngle = on; fArcFrac = frac; }
   /// A track is `contained` (stops in the gas) when its far-end centroid has radius < rMax and zMin < z_lab < zMax.
   /// The a1975 defaults: rMax 145 mm, z 10-955 mm -- tracks leaving the drift window end at z_lab ~964 mm otherwise
   /// look contained. Re-measure for another geometry / drift window.
   void SetContainment(double rMax, double zMin, double zMax) { fContRMax = rMax; fContZMin = zMin; fContZMax = zMax; }
   /// midPeak: charge profile (10 mm bins along the backbone from the vertex end, 3-bin smoothed) peaks before
   /// fpkMax of the track, more than minAfter mm before its end, and the mean after peak+20 mm is < tailMax x peak.
   void SetMidPeak(double fpkMax, double minAfter, double tailMax) { fMpFpk = fpkMax; fMpAfter = minAfter; fMpTail = tailMax; }

private:
   double fBField{2.85};
   double fZSign{-1}, fZOffset{1000};
   int fExpSense{1};
   int fMinHits{30};
   double fSmallPadRadius{152.0};
   double fBinW{4.0};
   double fMaxEdge{-4.0};
   int fMinBlob{8};
   double fBlobFrac{0.15}, fMaxDropFrac{0.5};
   double fBlobOffCircle{6.0}, fBlobLinMin{0.0};
   double fGuardFrac{0.6};
   double fDcaMax{10.0}, fRStartAxis{45.0};
   double fWinMaxFrac{1.0}, fWinMaxLen{1e9};
   bool fHitRefit{false}, fHitRefitZ{false};
   double fBraggThr{0.2}, fRContained{150.0};
   double fDedxGain{1.0};
   bool fArcAngle{false};
   double fArcFrac{0.5};
   double fContRMax{145.0}, fContZMin{10.0}, fContZMax{955.0};
   double fMpFpk{0.75}, fMpAfter{30.0}, fMpTail{0.35};
};

} // namespace AtTools
#endif
