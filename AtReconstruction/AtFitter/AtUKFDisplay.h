/// @file AtUKFDisplay.h
/// @brief Standalone interactive display for UKF track fitting.
///
/// Inspired by GenFit's EventDisplay: a self-contained singleton with
/// live refitting, 3D track visualization, and parameter controls.
/// Independent of AtViewerManager — uses TCanvas (no OpenGL/Eve required).
///
/// Usage:
///   auto *display = AtUKFDisplay::GetInstance();
///   display->LoadFiles("output_digi.root", "output_ukf_only.root");
///   display->GotoEvent(1);
///
/// Or from a macro:
///   root -l run_display.C

#ifndef ATUKFDISPLAY_H
#define ATUKFDISPLAY_H

#include <Rtypes.h>
#include <TNamed.h>

#include <memory>
#include <string>
#include <vector>

class TCanvas;
class TFile;
class TGCheckButton;
class TGMainFrame;
class TGNumberEntry;
class TGTextButton;
class TGComboBox;
class TGLabel;
class TPaveText;
class TTree;
class TClonesArray;
class TGraph;
class TGraphErrors;

class AtPatternEvent;
class AtTrackingEvent;
class AtFittedTrack;
class AtTrack;

namespace EventFit {
class AtFitterUKF;
}
namespace AtTools {
class AtELossModel;
}

class AtUKFDisplay : public TNamed {
public:
   /// Get the singleton instance. Creates it on first call.
   static AtUKFDisplay *GetInstance();

   /// Load input files.
   /// @param digiFile Path to output_digi.root (required: has AtPatternEvent)
   /// @param fittedFile Path to output_ukf_only.root (optional: has AtTrackingEvent)
   /// @param mcFile Path to attpcsim.root (optional: has MC truth for seeding)
   void LoadFiles(const char *digiFile, const char *fittedFile = nullptr, const char *mcFile = nullptr);

   /// Event navigation
   void GotoEvent(int id);
   void NextEvent(int step = 1);
   void PrevEvent(int step = 1);
   int GetNEvents() const { return fNEvents; }
   int GetCurrentEvent() const { return fCurrentEvent; }

   /// Fit the current track with current parameters
   /// @param skipReload If true, don't re-read tree (used after re-clustering)
   void FitCurrentTrack(bool skipReload = false);

   /// Display options: string of flags
   /// C=clusters, S=smoothed, F=filtered, R=residuals, E=errors
   void SetOptions(const char *opts) { fDrawOptions = opts; }

private:
   AtUKFDisplay();
   ~AtUKFDisplay();

   // GUI construction
   void MakeGui();
   void MakeControlPanel(TGMainFrame *mainFrame);
   void MakeDiagnosticsCanvas();

   // Drawing
   void DrawEvent();
   void DrawClusters(AtTrack &track);
   void DrawFittedTrack(const AtFittedTrack &fitted);
   void ClearEveElements();
   void UpdateDiagnostics(AtTrack &track, const AtFittedTrack *fitted);

   // Fitter management
   void CreateFitter();

public:
   // GUI callbacks (must be public for ROOT signal/slot)
   void GuiGotoEvent(); // *SIGNAL*
   void GuiNextEvent(); // *SIGNAL*
   void GuiPrevEvent(); // *SIGNAL*
   void GuiFit();       // *SIGNAL*
   void GuiRecluster(); // *SIGNAL*

private:

   // Singleton
   static AtUKFDisplay *fInstance;

   // Files and trees
   std::unique_ptr<TFile> fFileDigi;
   std::unique_ptr<TFile> fFileFit;
   std::unique_ptr<TFile> fFileMC;
   TTree *fTreeDigi{nullptr};
   TTree *fTreeFit{nullptr};
   TTree *fTreeMC{nullptr};
   TClonesArray *fPatEvtArr{nullptr};
   TClonesArray *fTrackingEvtArr{nullptr};
   TClonesArray *fMCTrackArr{nullptr};

   // Event state
   int fCurrentEvent{0};
   int fNEvents{0};
   double fZPadPlane{1000.0};

   // Display options
   std::string fDrawOptions{"CS"};
   bool fUseEve{false}; // True if TEveManager is available

   // GUI widgets
   TGMainFrame *fMainFrame{nullptr};
   TGNumberEntry *fEventEntry{nullptr};
   TGLabel *fEventInfoLabel{nullptr};
   TGNumberEntry *fAlphaEntry{nullptr};
   TGNumberEntry *fMeasSigmaEntry{nullptr};
   TGNumberEntry *fMinClustersEntry{nullptr};
   TGNumberEntry *fELossScaleEntry{nullptr};
   TGNumberEntry *fZPadPlaneEntry{nullptr};
   TGCheckButton *fStragglingBtn{nullptr};
   TGCheckButton *fPerClusterCovBtn{nullptr};
   TGCheckButton *fAutoFitBtn{nullptr};
   TGComboBox *fParticleBox{nullptr};
   TGNumberEntry *fMomSeedEntry{nullptr};
   TGNumberEntry *fMomSigmaEntry{nullptr};
   TGCheckButton *fUseMCTruthBtn{nullptr};
   TGCheckButton *fShowRawHitsBtn{nullptr};
   TGNumberEntry *fClusterRadiusEntry{nullptr};
   TGNumberEntry *fClusterDistEntry{nullptr};
   TGNumberEntry *fHitsPerClusterEntry{nullptr};
   TGComboBox *fClusterMethodBox{nullptr};

   // Canvases
   TCanvas *fTrackCanvas{nullptr}; // 2D projections (TCanvas mode)
   TCanvas *fDiagCanvas{nullptr};  // Diagnostics

   // Fitter (created on demand)
   std::unique_ptr<EventFit::AtFitterUKF> fFitter;

   ClassDef(AtUKFDisplay, 0);
};

#endif
