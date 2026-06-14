#ifndef ATTABFITTED_H
#define ATTABFITTED_H
#include "AtDataObserver.h"
#include "AtTabBase.h"
#include "AtViewerManagerSubject.h"

#include <Rtypes.h>
#include <TEveElement.h>
#include <TEveEventManager.h>
#include <TEveLine.h>
#include <TEvePointSet.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

class TBuffer;
class TClass;
class TClonesArray;
class TEveWindowSlot;
class TMemberInspector;
namespace DataHandling {
class AtSubject;
}

/**
 * @brief Tab that overlays AtFitterUKF results on the 3D view.
 *
 * For each AtFittedTrack draws the polyline (vertex, smoothed[0..N-1]) — the
 * trajectory points the UKF already integrated. No re-propagation, no helix
 * reconstruction; mirrors the approach in 16C_pp/AtUKFDisplay::DrawEvent.
 * Self-contained: looks up the AtTrackingEvent branch directly from
 * FairRootManager (bypasses the AtBranch index round-trip that resolved to
 * the wrong branch in chained-source setups).
 *
 * Display is in cm (matches AtTabMain's mm/10 hit rendering); SetFlipZ(true)
 * if the file's AtFittedTrack stores positions in lab z but the rest of
 * the display works in reco frame z = -lab.
 */
class AtTabFitted : public AtTabBase, public DataHandling::AtObserver {
protected:
   using TEvePointSetPtr = std::unique_ptr<TEvePointSet>;
   using TEveEventManagerPtr = std::unique_ptr<TEveEventManager>;

   TEveEventManagerPtr fEveTracking{std::make_unique<TEveEventManager>("AtTrackingEvent")};
   // TEveLine objects are owned by Eve once added; RemoveElements deletes
   // them, so the raw pointers we keep are dangling after a clear and we
   // never dereference them.
   std::vector<TEveLine *> fFittedTrackLines;
   // fVertexSet is persistent: pinned with SetDestroyOnZeroRefCnt(false) so
   // RemoveElements does not destroy it, and we re-AddElement each update.
   TEvePointSetPtr fVertexSet{std::make_unique<TEvePointSet>("Fitted vertex (Xtr)")};

   DataHandling::AtTreeEntry *fEntry{nullptr};
   TClonesArray *fTrackingArray{nullptr};

   double fZSign{1.0};
   bool fDrawVertex{true};
   bool fInitialized{false};

   std::string fBranchName{"AtTrackingEvent"}; //< Branch to read fitted tracks from
   bool fDrawSmoothed{false};                  //< Draw smoothed polyline instead of analytical helix
   Color_t fFixedColor{-1};                    //< <0: per-track palette; else single colour
   Style_t fLineStyle{1};                      //< Line style (1=solid, 2=dashed, ...)
   int fLineWidthCfg{2};                       //< Line width
   double fZPadPlane{0.0};                      //< Pad-plane z (mm) for lab->digi z mapping (smoothed mode)

   // Analytical helix sampling (from fit kinematics)
   double fBField_T{4.0};      //< Solenoid B (Tesla); along +ẑ for PUMA.
   double fMaxArc_mm{500.0};   //< Stop the helix after this arc length.
   int fHelixSamples{200};     //< Number of points sampled along the arc.

public:
   explicit AtTabFitted(const std::string &name = "Fitted", const std::string &branchName = "AtTrackingEvent");
   ~AtTabFitted() override;
   void InitTab() override;
   void Exec() override {}
   void Update(DataHandling::AtSubject *sub) override;

   /// Flip z when drawing fitted trajectory points. Set true when
   /// AtFittedTrack stores lab z but the display is in reco frame.
   void SetFlipZ(bool flip) { fZSign = flip ? -1.0 : 1.0; }
   void SetDrawVertex(bool draw) { fDrawVertex = draw; }
   /// Solenoid B-field along z (Tesla). Default 4 T (PUMA).
   void SetBField(double Bz_T) { fBField_T = Bz_T; }
   /// Stop the propagated helix after this arc length (mm) and use this
   /// many sample points along it. 500 mm / 200 samples is enough to draw
   /// a smooth curve through the PUMA TPC volume.
   void SetHelixSampling(double maxArc_mm, int nSamples)
   {
      fMaxArc_mm = maxArc_mm;
      fHelixSamples = nSamples;
   }

   /// Read fitted tracks from this branch (default "AtTrackingEvent"). Point two
   /// AtTabFitted tabs at two fitters' outputs (e.g. "AtTrackingEventUKF" and
   /// "AtTrackingEventGenfit") to overlay both Kalman filters in one view.
   void SetBranchName(const std::string &name) { fBranchName = name; }
   /// Draw the actual fitted trajectory (vertex + GetSmoothedPositions polyline)
   /// instead of an analytical helix reconstructed from the fit kinematics. This is
   /// the honest per-filter trajectory and is fitter-agnostic (both UKF and genfit
   /// store it via SetSmoothedPositions). Default false (helix, legacy PUMA behaviour).
   void SetDrawSmoothed(bool s) { fDrawSmoothed = s; }
   /// Force a single colour for every track of this tab (e.g. blue=UKF, red=genfit).
   /// <0 keeps the per-track palette.
   void SetTrackColor(Color_t c) { fFixedColor = c; }
   void SetLineStyle(Style_t s) { fLineStyle = s; }
   void SetLineWidth(int w) { fLineWidthCfg = w; }
   /// Pad-plane z (mm). When >0 and drawing smoothed trajectories, the fitted
   /// lab-frame z (= ZPadPlane - z_digi as stored by the fitters) is mapped back to
   /// the digi frame z_digi = ZPadPlane - z_lab so the trajectory overlays the
   /// (digi-frame) hits drawn by AtTabMain. 0 keeps the SetFlipZ behaviour.
   void SetZPadPlane(double z) { fZPadPlane = z; }

protected:
   void MakeTab(TEveWindowSlot * /*slot*/) override {}

   void UpdateFittedElements();

   ClassDefOverride(AtTabFitted, 1)
};

#endif
