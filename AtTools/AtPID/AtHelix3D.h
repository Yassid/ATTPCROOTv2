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
   void SetDedxGain(double g) { fDedxGain = g; }
   double GetDedxGain() const { return fDedxGain; }

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
   double fBraggThr{0.2}, fRContained{150.0};
   double fDedxGain{1.0};
};

} // namespace AtTools
#endif
