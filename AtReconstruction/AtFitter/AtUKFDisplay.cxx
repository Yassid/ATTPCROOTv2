#include "AtUKFDisplay.h"

#include "AtELossCATIMA.h"
#include "AtFittedTrack.h"
#include "AtFitterUKF.h"
#include "AtHitCluster.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtMCTrack.h"
#include "AtTrackTransformer.h"
#include "AtTrackingEvent.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGComboBox.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGLayout.h>
#include <TGNumberEntry.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <memory>
#include <tuple>

ClassImp(AtUKFDisplay);

AtUKFDisplay *AtUKFDisplay::fInstance = nullptr;

// ===========================================================================
// Singleton
// ===========================================================================
AtUKFDisplay *AtUKFDisplay::GetInstance()
{
   if (!fInstance)
      fInstance = new AtUKFDisplay();
   return fInstance;
}

AtUKFDisplay::AtUKFDisplay() : TNamed("AtUKFDisplay", "UKF Fitter Display")
{
   std::cout << "AtUKFDisplay: creating..." << std::endl;
   fUseEve = false; // TCanvas mode only — no Eve dependency

   MakeGui();
   MakeDiagnosticsCanvas();
   std::cout << "AtUKFDisplay: ready." << std::endl;
}

AtUKFDisplay::~AtUKFDisplay()
{
   fInstance = nullptr;
}

// ===========================================================================
// File loading
// ===========================================================================
void AtUKFDisplay::LoadFiles(const char *digiFile, const char *fittedFile, const char *mcFile)
{
   fFileDigi.reset(TFile::Open(digiFile));
   if (!fFileDigi || fFileDigi->IsZombie()) {
      std::cerr << "AtUKFDisplay: cannot open " << digiFile << std::endl;
      return;
   }
   fTreeDigi = (TTree *)fFileDigi->Get("cbmsim");
   if (!fTreeDigi) {
      std::cerr << "AtUKFDisplay: no cbmsim tree in " << digiFile << std::endl;
      return;
   }
   fTreeDigi->SetBranchAddress("AtPatternEvent", &fPatEvtArr);
   fNEvents = fTreeDigi->GetEntries();
   std::cout << "AtUKFDisplay: loaded " << fNEvents << " events from " << digiFile << std::endl;

   if (fittedFile) {
      fFileFit.reset(TFile::Open(fittedFile));
      if (fFileFit && !fFileFit->IsZombie()) {
         fTreeFit = (TTree *)fFileFit->Get("cbmsim");
         if (fTreeFit) {
            fTreeFit->SetBranchAddress("AtTrackingEvent", &fTrackingEvtArr);
            std::cout << "AtUKFDisplay: fitted data loaded from " << fittedFile << std::endl;
         } else {
            std::cerr << "AtUKFDisplay: no cbmsim tree in " << fittedFile << " (skipping)" << std::endl;
            fTreeFit = nullptr;
         }
      } else {
         std::cerr << "AtUKFDisplay: cannot open " << fittedFile << " (skipping, use Fit button)" << std::endl;
         fFileFit.reset();
         fTreeFit = nullptr;
      }
   }

   if (mcFile) {
      fFileMC.reset(TFile::Open(mcFile));
      if (fFileMC && !fFileMC->IsZombie()) {
         fTreeMC = (TTree *)fFileMC->Get("cbmsim");
         if (fTreeMC) {
            fTreeMC->SetBranchAddress("MCTrack", &fMCTrackArr);
            std::cout << "AtUKFDisplay: MC truth loaded from " << mcFile << std::endl;
         }
      }
   }
}

// ===========================================================================
// Event navigation
// ===========================================================================
void AtUKFDisplay::GotoEvent(int id)
{
   if (id < 0 || id >= fNEvents)
      return;
   fCurrentEvent = id;
   if (fEventEntry)
      fEventEntry->SetIntNumber(id);
   DrawEvent();
}

void AtUKFDisplay::NextEvent(int step)
{
   // Find next event with a track
   for (int i = fCurrentEvent + step; i < fNEvents; i++) {
      fTreeDigi->GetEntry(i);
      if (fPatEvtArr && fPatEvtArr->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)fPatEvtArr->At(0);
         if (!pe->GetTrackCand().empty() && pe->GetTrackCand()[0].GetHitClusterArray()->size() > 10) {
            GotoEvent(i);
            return;
         }
      }
   }
   std::cout << "No more events with tracks." << std::endl;
}

void AtUKFDisplay::PrevEvent(int step)
{
   for (int i = fCurrentEvent - step; i >= 0; i--) {
      fTreeDigi->GetEntry(i);
      if (fPatEvtArr && fPatEvtArr->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)fPatEvtArr->At(0);
         if (!pe->GetTrackCand().empty() && pe->GetTrackCand()[0].GetHitClusterArray()->size() > 10) {
            GotoEvent(i);
            return;
         }
      }
   }
   std::cout << "No previous events with tracks." << std::endl;
}

// ===========================================================================
// GUI construction
// ===========================================================================
void AtUKFDisplay::MakeGui()
{
   fMainFrame = new TGMainFrame(gClient->GetRoot(), 280, 700);
   fMainFrame->SetWindowName("UKF Fitter Controls");
   fMainFrame->SetCleanup(kDeepCleanup);

   MakeControlPanel(fMainFrame);

   fMainFrame->MapSubwindows();
   fMainFrame->Resize(280, 700);
   fMainFrame->MapWindow();
}

void AtUKFDisplay::MakeControlPanel(TGMainFrame *mf)
{
   auto *vf = new TGVerticalFrame(mf);

   // --- Event navigation ---
   auto *evtFrame = new TGHorizontalFrame(vf);
   auto *prevBtn = new TGTextButton(evtFrame, " << ");
   prevBtn->Connect("Clicked()", "AtUKFDisplay", this, "GuiPrevEvent()");
   evtFrame->AddFrame(prevBtn, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

   fEventEntry = new TGNumberEntry(evtFrame, 0, 6, -1, TGNumberFormat::kNESInteger);
   evtFrame->AddFrame(fEventEntry, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   auto *nextBtn = new TGTextButton(evtFrame, " >> ");
   nextBtn->Connect("Clicked()", "AtUKFDisplay", this, "GuiNextEvent()");
   evtFrame->AddFrame(nextBtn, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));

   vf->AddFrame(evtFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 2));

   auto *goBtn = new TGTextButton(vf, "Go to Event");
   goBtn->Connect("Clicked()", "AtUKFDisplay", this, "GuiGotoEvent()");
   vf->AddFrame(goBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   // --- Event info ---
   fEventInfoLabel = new TGLabel(vf, "No event loaded");
   vf->AddFrame(fEventInfoLabel, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 2));

   // --- Separator ---
   vf->AddFrame(new TGHorizontalFrame(vf, 1, 2), new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 5));

   // --- Fitter parameters ---
   auto *parLabel = new TGLabel(vf, "=== Fitter Parameters ===");
   vf->AddFrame(parLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

   // Particle selection
   auto *partFrame = new TGHorizontalFrame(vf);
   partFrame->AddFrame(new TGLabel(partFrame, "Particle:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fParticleBox = new TGComboBox(partFrame, -1);
   fParticleBox->AddEntry("Proton", 0);
   fParticleBox->AddEntry("Deuteron", 1);
   fParticleBox->AddEntry("Alpha", 2);
   fParticleBox->AddEntry("Pi+", 3);
   fParticleBox->AddEntry("Pi-", 4);
   fParticleBox->Select(0);
   fParticleBox->Resize(100, 20);
   partFrame->AddFrame(fParticleBox, new TGLayoutHints(kLHintsRight | kLHintsExpandX, 2, 2, 2, 2));
   vf->AddFrame(partFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   // Alpha
   auto *alphaFrame = new TGHorizontalFrame(vf);
   alphaFrame->AddFrame(new TGLabel(alphaFrame, "Alpha:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fAlphaEntry = new TGNumberEntry(alphaFrame, 1e-3, 8, -1, TGNumberFormat::kNESRealFour);
   alphaFrame->AddFrame(fAlphaEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(alphaFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Measurement sigma
   auto *sigFrame = new TGHorizontalFrame(vf);
   sigFrame->AddFrame(new TGLabel(sigFrame, "Meas sigma [mm]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fMeasSigmaEntry = new TGNumberEntry(sigFrame, 2.0, 6, -1, TGNumberFormat::kNESRealTwo);
   sigFrame->AddFrame(fMeasSigmaEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(sigFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Min clusters
   auto *minFrame = new TGHorizontalFrame(vf);
   minFrame->AddFrame(new TGLabel(minFrame, "Min clusters:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fMinClustersEntry = new TGNumberEntry(minFrame, 10, 4, -1, TGNumberFormat::kNESInteger);
   minFrame->AddFrame(fMinClustersEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(minFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // E-loss scale
   auto *elossFrame = new TGHorizontalFrame(vf);
   elossFrame->AddFrame(new TGLabel(elossFrame, "E-loss scale:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fELossScaleEntry = new TGNumberEntry(elossFrame, 1.0, 6, -1, TGNumberFormat::kNESRealTwo);
   elossFrame->AddFrame(fELossScaleEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(elossFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // ZPadPlane
   auto *zpFrame = new TGHorizontalFrame(vf);
   zpFrame->AddFrame(new TGLabel(zpFrame, "ZPadPlane [mm]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fZPadPlaneEntry = new TGNumberEntry(zpFrame, 1000.0, 6, -1, TGNumberFormat::kNESRealOne);
   zpFrame->AddFrame(fZPadPlaneEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(zpFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Momentum seed override (0 = use Brho)
   auto *momFrame = new TGHorizontalFrame(vf);
   momFrame->AddFrame(new TGLabel(momFrame, "p seed [MeV/c]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fMomSeedEntry = new TGNumberEntry(momFrame, 0, 6, -1, TGNumberFormat::kNESRealOne);
   momFrame->AddFrame(fMomSeedEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(momFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Momentum sigma fraction (initial covariance)
   auto *sigFracFrame = new TGHorizontalFrame(vf);
   sigFracFrame->AddFrame(new TGLabel(sigFracFrame, "p sigma frac:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fMomSigmaEntry = new TGNumberEntry(sigFracFrame, 0.3, 6, -1, TGNumberFormat::kNESRealTwo);
   sigFracFrame->AddFrame(fMomSigmaEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(sigFracFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // B-field Bz [T]
   auto *bFrame = new TGHorizontalFrame(vf);
   bFrame->AddFrame(new TGLabel(bFrame, "Bz [T]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fBFieldEntry = new TGNumberEntry(bFrame, 2.85, 6, -1, TGNumberFormat::kNESRealTwo);
   bFrame->AddFrame(fBFieldEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(bFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Gas density [g/cm^3] (CATIMA)
   auto *rhoFrame = new TGHorizontalFrame(vf);
   rhoFrame->AddFrame(new TGLabel(rhoFrame, "Gas rho [g/cm3]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fGasDensityEntry = new TGNumberEntry(rhoFrame, 3.553e-5, 10, -1, TGNumberFormat::kNESRealFour);
   rhoFrame->AddFrame(fGasDensityEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(rhoFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // === Clustering controls ===
   vf->AddFrame(new TGHorizontalFrame(vf, 1, 2), new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 5));
   auto *clLabel = new TGLabel(vf, "=== Clustering ===");
   vf->AddFrame(clLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

   // Clustering method
   auto *methFrame = new TGHorizontalFrame(vf);
   methFrame->AddFrame(new TGLabel(methFrame, "Method:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fClusterMethodBox = new TGComboBox(methFrame, -1);
   fClusterMethodBox->AddEntry("Smooth3D", 0);
   fClusterMethodBox->AddEntry("GroupByN", 1);
   fClusterMethodBox->Select(0);
   fClusterMethodBox->Resize(100, 20);
   methFrame->AddFrame(fClusterMethodBox, new TGLayoutHints(kLHintsRight | kLHintsExpandX, 2, 2, 2, 2));
   vf->AddFrame(methFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // Smooth3D params
   auto *radFrame = new TGHorizontalFrame(vf);
   radFrame->AddFrame(new TGLabel(radFrame, "Radius [mm]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fClusterRadiusEntry = new TGNumberEntry(radFrame, 10.0, 6, -1, TGNumberFormat::kNESRealOne);
   radFrame->AddFrame(fClusterRadiusEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(radFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   auto *distFrame = new TGHorizontalFrame(vf);
   distFrame->AddFrame(new TGLabel(distFrame, "Distance [mm]:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fClusterDistEntry = new TGNumberEntry(distFrame, 20.0, 6, -1, TGNumberFormat::kNESRealOne);
   distFrame->AddFrame(fClusterDistEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(distFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   // GroupByN param
   auto *hpcFrame = new TGHorizontalFrame(vf);
   hpcFrame->AddFrame(new TGLabel(hpcFrame, "Hits/cluster:"), new TGLayoutHints(kLHintsLeft, 2, 2, 3, 2));
   fHitsPerClusterEntry = new TGNumberEntry(hpcFrame, 15, 4, -1, TGNumberFormat::kNESInteger);
   hpcFrame->AddFrame(fHitsPerClusterEntry, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   vf->AddFrame(hpcFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));

   fShowRawHitsBtn = new TGCheckButton(vf, "Show raw hits");
   vf->AddFrame(fShowRawHitsBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 3, 1));

   auto *reclusterBtn = new TGTextButton(vf, "  Re-cluster + Fit  ");
   reclusterBtn->Connect("Clicked()", "AtUKFDisplay", this, "GuiRecluster()");
   reclusterBtn->SetBackgroundColor(0x0088ff);
   vf->AddFrame(reclusterBtn, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 3));

   // === Fitter options ===
   vf->AddFrame(new TGHorizontalFrame(vf, 1, 2), new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 5));

   // Checkboxes
   fUseMCTruthBtn = new TGCheckButton(vf, "Use MC truth seed");
   vf->AddFrame(fUseMCTruthBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 3, 1));

   fStragglingBtn = new TGCheckButton(vf, "Energy straggling");
   vf->AddFrame(fStragglingBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 1, 1));

   fPerClusterCovBtn = new TGCheckButton(vf, "Per-cluster covariance");
   vf->AddFrame(fPerClusterCovBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 1, 1));

   fAutoFitBtn = new TGCheckButton(vf, "Auto-fit on navigate");
   vf->AddFrame(fAutoFitBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 1, 1));

   // --- Separator ---
   vf->AddFrame(new TGHorizontalFrame(vf, 1, 2), new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 5));

   // Fit button
   auto *fitBtn = new TGTextButton(vf, "  Fit Current Track  ");
   fitBtn->Connect("Clicked()", "AtUKFDisplay", this, "GuiFit()");
   fitBtn->SetBackgroundColor(0x00cc00); // Green
   vf->AddFrame(fitBtn, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   mf->AddFrame(vf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
}

void AtUKFDisplay::MakeDiagnosticsCanvas()
{
   fTrackCanvas = new TCanvas("cUKFTrack", "UKF Track Display", 1200, 500);
   fTrackCanvas->Divide(3, 1);

   fDiagCanvas = new TCanvas("cUKFDiag", "UKF Diagnostics", 1000, 500);
   fDiagCanvas->Divide(4, 1);
}

// ===========================================================================
// Drawing
// ===========================================================================
void AtUKFDisplay::ClearEveElements()
{
   // No-op in TCanvas mode
}

void AtUKFDisplay::DrawEvent()
{
   if (!fTreeDigi)
      return;

   fTreeDigi->GetEntry(fCurrentEvent);
   if (fTreeFit)
      fTreeFit->GetEntry(fCurrentEvent);

   ClearEveElements();

   if (!fPatEvtArr || fPatEvtArr->GetEntries() == 0) {
      if (fEventInfoLabel)
         fEventInfoLabel->SetText(Form("Event %d: no data", fCurrentEvent));
      // No Eve redraw needed in TCanvas mode
      return;
   }

   auto *patEvt = (AtPatternEvent *)fPatEvtArr->At(0);
   auto &tracks = patEvt->GetTrackCand();
   if (tracks.empty()) {
      if (fEventInfoLabel)
         fEventInfoLabel->SetText(Form("Event %d: no tracks", fCurrentEvent));
      // No Eve redraw needed in TCanvas mode
      return;
   }

   // Use the largest track
   int bestTrack = 0;
   for (size_t t = 1; t < tracks.size(); t++) {
      if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
         bestTrack = t;
   }
   auto &track = tracks[bestTrack];

   // Get fitted track if available
   AtFittedTrack *fitted = nullptr;
   if (fTrackingEvtArr && fTrackingEvtArr->GetEntries() > 0) {
      auto *te = (AtTrackingEvent *)fTrackingEvtArr->At(0);
      auto &ft = te->GetFittedTracks();
      if (!ft.empty())
         fitted = ft[0].get();
   }

   int nClusters = track.GetHitClusterArray()->size();

   if (fTrackCanvas) {
      // Raw hits (before clustering)
      std::vector<double> hx, hy, hz;
      bool showRaw = fShowRawHitsBtn && fShowRawHitsBtn->IsOn();
      if (showRaw) {
         auto &hitArr = track.GetHitArray();
         for (auto &hit : hitArr) {
            auto pos = hit->GetPosition();
            hx.push_back(pos.X());
            hy.push_back(pos.Y());
            hz.push_back(fZPadPlane - pos.Z());
         }
      }

      // Clusters
      auto *clusters = track.GetHitClusterArray();
      std::vector<double> cx, cy, cz;
      for (auto &cl : *clusters) {
         auto pos = cl.GetPosition();
         cx.push_back(pos.X());
         cy.push_back(pos.Y());
         cz.push_back(fZPadPlane - pos.Z());
      }

      // Smoothed positions
      std::vector<double> sx, sy, sz;
      if (fitted) {
         auto vtx = fitted->GetVertex();
         sx.push_back(vtx.X());
         sy.push_back(vtx.Y());
         sz.push_back(vtx.Z());
         for (auto &sp : fitted->GetSmoothedPositions()) {
            sx.push_back(sp.X());
            sy.push_back(sp.Y());
            sz.push_back(sp.Z());
         }
      }

      // Helper lambda to draw one projection pad with raw hits, clusters, and fit
      auto drawPad = [&](int padNum, const char *title,
                         std::vector<double> &ax, std::vector<double> &ay,   // clusters
                         std::vector<double> &bx, std::vector<double> &by,   // raw hits
                         std::vector<double> &fx, std::vector<double> &fy) { // fitted
         fTrackCanvas->cd(padNum);
         gPad->Clear();

         // Draw raw hits first (background)
         TGraph *gRaw = nullptr;
         if (showRaw && !bx.empty()) {
            gRaw = new TGraph(bx.size(), bx.data(), by.data());
            gRaw->SetTitle(title);
            gRaw->SetMarkerStyle(1);
            gRaw->SetMarkerSize(1);
            gRaw->SetMarkerColor(kGray + 1);
            gRaw->Draw("AP");
         }

         // Draw clusters
         auto *gCl = new TGraph(ax.size(), ax.data(), ay.data());
         if (!gRaw)
            gCl->SetTitle(title);
         gCl->SetMarkerStyle(20);
         gCl->SetMarkerSize(0.6);
         gCl->SetMarkerColor(kBlue);
         gCl->Draw(gRaw ? "P SAME" : "AP");

         // Draw fitted
         TGraph *gFit = nullptr;
         if (!fx.empty()) {
            gFit = new TGraph(fx.size(), fx.data(), fy.data());
            gFit->SetMarkerStyle(24);
            gFit->SetMarkerSize(0.5);
            gFit->SetMarkerColor(kRed);
            gFit->Draw("P SAME");
         }

         // Legend
         auto *leg = new TLegend(0.12, 0.70, 0.50, 0.88);
         leg->SetTextSize(0.03);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         if (gRaw)
            leg->AddEntry(gRaw, Form("Raw hits (%d)", (int)bx.size()), "p");
         leg->AddEntry(gCl, Form("Clusters (%d)", (int)ax.size()), "p");
         if (gFit)
            leg->AddEntry(gFit, "Fitted (UKF)", "p");
         leg->Draw();
      };

      drawPad(1, Form("Event %d  XY;X [mm];Y [mm]", fCurrentEvent), cx, cy, hx, hy, sx, sy);
      drawPad(2, "XZ;Z_{lab} [mm];X [mm]", cz, cx, hz, hx, sz, sx);
      drawPad(3, "YZ;Z_{lab} [mm];Y [mm]", cz, cy, hz, hy, sz, sy);

      fTrackCanvas->Update();
   }

   // Update info label
   if (fitted) {
      auto kin = fitted->GetKinematics();
      double p = std::sqrt(2 * 938.272 * kin.kineticEnergy + kin.kineticEnergy * kin.kineticEnergy);
      if (fEventInfoLabel)
         fEventInfoLabel->SetText(
            Form("Ev %d: %d cl, p=%.1f MeV/c, KE=%.2f MeV", fCurrentEvent, nClusters, p, kin.kineticEnergy));
   } else {
      if (fEventInfoLabel)
         fEventInfoLabel->SetText(Form("Ev %d: %d clusters, no fit", fCurrentEvent, nClusters));
   }

   // Update diagnostics
   UpdateDiagnostics(track, fitted);

   // Auto-fit if enabled
   if (fAutoFitBtn && fAutoFitBtn->IsOn() && !fitted)
      FitCurrentTrack();
}

void AtUKFDisplay::DrawClusters(AtTrack & /*track*/)
{
   // Clusters are drawn in DrawEvent via TCanvas
}

void AtUKFDisplay::DrawFittedTrack(const AtFittedTrack & /*fitted*/)
{
   // Fitted track is drawn in DrawEvent via TCanvas
}

void AtUKFDisplay::UpdateDiagnostics(AtTrack &track, const AtFittedTrack *fitted)
{
   if (!fDiagCanvas)
      return;

   auto *clusters = track.GetHitClusterArray();
   int nCl = clusters->size();

   // Cluster charge profile
   fDiagCanvas->cd(1);
   gPad->Clear();
   auto *gQ = new TGraph(nCl);
   for (int i = 0; i < nCl; i++)
      gQ->SetPoint(i, i, clusters->at(i).GetCharge());
   gQ->SetTitle("Charge;Cluster;Q [arb]");
   gQ->SetMarkerStyle(20);
   gQ->SetMarkerSize(0.4);
   gQ->SetMarkerColor(kGreen + 2);
   gQ->Draw("AP");

   if (fitted) {
      auto &smoothed = fitted->GetSmoothedPositions();
      auto vtx = fitted->GetVertex();
      auto kin = fitted->GetKinematics();

      // Build smoothed position list
      std::vector<ROOT::Math::XYZPoint> sPts;
      sPts.emplace_back(vtx.X(), vtx.Y(), vtx.Z());
      for (auto &sp : smoothed)
         sPts.push_back(sp);

      // Residuals
      fDiagCanvas->cd(2);
      gPad->Clear();
      int nPts = std::min((int)sPts.size(), nCl);
      auto *gResid = new TGraph(nPts);
      for (int i = 0; i < nPts; i++) {
         ROOT::Math::XYZPoint cp(clusters->at(i).GetPosition().X(), clusters->at(i).GetPosition().Y(),
                                 fZPadPlane - clusters->at(i).GetPosition().Z());
         double d = (sPts[i] - cp).R();
         gResid->SetPoint(i, i, d);
      }
      gResid->SetTitle("Residuals;Cluster;|sm-meas| [mm]");
      gResid->SetMarkerStyle(20);
      gResid->SetMarkerSize(0.4);
      gResid->SetMarkerColor(kMagenta);
      gResid->Draw("AP");

      // Z profile (arc length)
      fDiagCanvas->cd(3);
      gPad->Clear();
      auto *gZ = new TGraph(nCl);
      for (int i = 0; i < nCl; i++)
         gZ->SetPoint(i, i, fZPadPlane - clusters->at(i).GetPosition().Z());
      gZ->SetTitle("Z_{lab};Cluster;Z [mm]");
      gZ->SetMarkerStyle(20);
      gZ->SetMarkerSize(0.4);
      gZ->Draw("AP");

      // Fit info
      fDiagCanvas->cd(4);
      gPad->Clear();
      double p = std::sqrt(2 * 938.272 * kin.kineticEnergy + kin.kineticEnergy * kin.kineticEnergy);
      auto *pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
      pt->SetTextAlign(12);
      pt->SetTextSize(0.07);
      pt->SetFillColor(0);
      pt->AddText(Form("Event %d", fCurrentEvent));
      pt->AddText(Form("Clusters: %d", nCl));
      pt->AddText(Form("p = %.1f MeV/c", p));
      pt->AddText(Form("KE = %.2f MeV", kin.kineticEnergy));
      pt->AddText(Form("#theta = %.1f#circ", kin.theta * 180.0 / TMath::Pi()));
      pt->AddText(Form("#phi = %.1f#circ", kin.phi * 180.0 / TMath::Pi()));
      pt->Draw();
   } else {
      fDiagCanvas->cd(4);
      gPad->Clear();
      auto *pt = new TPaveText(0.1, 0.3, 0.9, 0.7, "NDC");
      pt->SetTextSize(0.08);
      pt->AddText(Form("Event %d", fCurrentEvent));
      pt->AddText("No fit result");
      pt->Draw();
   }

   fDiagCanvas->Update();
}

// ===========================================================================
// Fitting
// ===========================================================================
void AtUKFDisplay::SetParticle(int kind)
{
   if (fParticleBox) fParticleBox->Select(kind);
}
void AtUKFDisplay::SetBField(double bz_T)
{
   if (fBFieldEntry) fBFieldEntry->SetNumber(bz_T);
}
void AtUKFDisplay::SetGasDensity(double rho_gcc)
{
   if (fGasDensityEntry) fGasDensityEntry->SetNumber(rho_gcc);
}

void AtUKFDisplay::CreateFitter()
{
   double charge = 1.602176634e-19;
   double mass = 938.272;

   // Check particle selection
   if (fParticleBox) {
      int sel = fParticleBox->GetSelected();
      if (sel == 1) { // Deuteron
         mass = 1875.613;
      } else if (sel == 2) { // Alpha
         charge *= 2;
         mass = 3727.379;
      } else if (sel == 3) { // Pi+
         mass = 139.57039;
      } else if (sel == 4) { // Pi-
         mass = 139.57039;
         charge = -charge;
      }
   }

   const double rho = fGasDensityEntry ? fGasDensityEntry->GetNumber() : 3.553e-5;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(rho);
   const double mass_amu = mass / 931.49401;
   eloss->SetProjectile(1, 1, mass_amu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   fFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass, std::move(eloss));
   const double bz = fBFieldEntry ? fBFieldEntry->GetNumber() : 2.85;
   fFitter->SetBField({0, 0, bz});

   // Read parameters from GUI
   double alpha = fAlphaEntry ? fAlphaEntry->GetNumber() : 1e-3;
   double measSigma = fMeasSigmaEntry ? fMeasSigmaEntry->GetNumber() : 2.0;
   int minClusters = fMinClustersEntry ? fMinClustersEntry->GetIntNumber() : 10;
   double eLossScale = fELossScaleEntry ? fELossScaleEntry->GetNumber() : 1.0;
   double zPadPlane = fZPadPlaneEntry ? fZPadPlaneEntry->GetNumber() : 1000.0;
   bool straggling = fStragglingBtn ? fStragglingBtn->IsOn() : false;
   bool perClusterCov = fPerClusterCovBtn ? fPerClusterCovBtn->IsOn() : false;

   fFitter->SetUKFParameters(alpha, 2.0, 0.0);
   fFitter->SetMeasurementSigma(measSigma);
   fFitter->SetMinClusters(minClusters);
   fFitter->SetEnableEnergyStraggling(straggling);
   fFitter->SetUsePerClusterCov(perClusterCov);
   fFitter->SetZPadPlane(zPadPlane);
   double momSigmaFrac = fMomSigmaEntry ? fMomSigmaEntry->GetNumber() : 0.3;
   fFitter->SetMomentumSigmaFrac(momSigmaFrac);

   // Momentum seed override (0 = use Brho from pattern recognition)
   double momSeed = fMomSeedEntry ? fMomSeedEntry->GetNumber() : 0;
   if (momSeed > 0)
      fFitter->SetMomentumSeed(momSeed);

   fZPadPlane = zPadPlane;
}

void AtUKFDisplay::FitCurrentTrack(bool skipReload)
{
   if (!fTreeDigi)
      return;

   if (!skipReload)
      fTreeDigi->GetEntry(fCurrentEvent);
   if (!fPatEvtArr || fPatEvtArr->GetEntries() == 0)
      return;

   auto *patEvt = (AtPatternEvent *)fPatEvtArr->At(0);
   auto &tracks = patEvt->GetTrackCand();
   if (tracks.empty())
      return;

   // Recreate fitter with current GUI parameters
   CreateFitter();

   // If MC truth seed is enabled, find the proton and override the seed
   bool useMCTruth = fUseMCTruthBtn && fUseMCTruthBtn->IsOn() && fTreeMC;
   if (useMCTruth) {
      fTreeMC->GetEntry(fCurrentEvent);
      if (fMCTrackArr) {
         for (int j = 0; j < fMCTrackArr->GetEntries(); j++) {
            auto *mcTrack = (AtMCTrack *)fMCTrackArr->At(j);
            if (mcTrack->GetPdgCode() == 2212) { // proton
               double pMC = std::sqrt(mcTrack->GetPx() * mcTrack->GetPx() + mcTrack->GetPy() * mcTrack->GetPy() +
                                      mcTrack->GetPz() * mcTrack->GetPz()) *
                            1e3; // GeV → MeV
               fFitter->SetMomentumSeed(pMC);
               std::cout << "  MC truth seed: p=" << pMC << " MeV/c" << std::endl;
               break;
            }
         }
      }
   }

   fFitter->Init();

   // Fit using the public FitEvent interface
   AtTrackingEvent trackingEvent;
   fFitter->FitEvent(&trackingEvent, patEvt);

   auto &fittedTracks = trackingEvent.GetFittedTracks();
   if (fittedTracks.empty()) {
      std::cout << "AtUKFDisplay: fit failed for event " << fCurrentEvent << std::endl;
      return;
   }

   auto kin = fittedTracks[0]->GetKinematics();
   double pFit = std::sqrt(2 * 938.272 * kin.kineticEnergy + kin.kineticEnergy * kin.kineticEnergy);

   // Print fit summary with initial covariance info
   double momSeed = fMomSeedEntry ? fMomSeedEntry->GetNumber() : 0;
   double momSigmaFrac = fMomSigmaEntry ? fMomSigmaEntry->GetNumber() : 0.3;
   double pSeed = momSeed > 0 ? momSeed : 37.0; // approximate Brho value
   std::cout << "=== UKF Fit Result ===" << std::endl;
   std::cout << "  Event:    " << fCurrentEvent << std::endl;
   std::cout << "  Seed:     p=" << pSeed << " MeV/c" << std::endl;
   std::cout << "  Init Cov: sigma_pos=" << 2.0 << " mm, sigma_p=" << momSigmaFrac * pSeed
             << " MeV/c (" << momSigmaFrac * 100 << "%), sigma_ang=" << 5.0 << " deg" << std::endl;
   std::cout << "  Result:   p=" << pFit << " MeV/c, KE=" << kin.kineticEnergy
             << " MeV, theta=" << kin.theta * 180.0 / TMath::Pi() << " deg" << std::endl;
   std::cout << "======================" << std::endl;

   // Use largest track for display
   int bestTrack = 0;
   for (size_t t = 1; t < tracks.size(); t++) {
      if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
         bestTrack = t;
   }

   // Redraw track canvas with fit overlaid
   if (fTrackCanvas) {
      auto *clusters = tracks[bestTrack].GetHitClusterArray();
      std::vector<double> cx, cy, cz;
      for (auto &cl : *clusters) {
         auto pos = cl.GetPosition();
         cx.push_back(pos.X());
         cy.push_back(pos.Y());
         cz.push_back(fZPadPlane - pos.Z());
      }

      std::vector<double> sx, sy, sz;
      auto vtx = fittedTracks[0]->GetVertex();
      sx.push_back(vtx.X());
      sy.push_back(vtx.Y());
      sz.push_back(vtx.Z());
      for (auto &sp : fittedTracks[0]->GetSmoothedPositions()) {
         sx.push_back(sp.X());
         sy.push_back(sp.Y());
         sz.push_back(sp.Z());
      }

      auto drawLegend = [](TGraph *gCl, TGraph *gFit) {
         auto *leg = new TLegend(0.12, 0.75, 0.45, 0.88);
         leg->SetTextSize(0.035);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         leg->AddEntry(gCl, "Clusters", "p");
         if (gFit)
            leg->AddEntry(gFit, "Fitted (UKF)", "p");
         leg->Draw();
      };

      fTrackCanvas->cd(1);
      gPad->Clear();
      auto *g1 = new TGraph(cx.size(), cx.data(), cy.data());
      g1->SetTitle(Form("Event %d XY (live fit);X [mm];Y [mm]", fCurrentEvent));
      g1->SetMarkerStyle(20);
      g1->SetMarkerSize(0.6);
      g1->SetMarkerColor(kBlue);
      g1->Draw("AP");
      auto *g1f = new TGraph(sx.size(), sx.data(), sy.data());
      g1f->SetMarkerStyle(24);
      g1f->SetMarkerSize(0.5);
      g1f->SetMarkerColor(kRed);
      g1f->Draw("P SAME");
      drawLegend(g1, g1f);

      fTrackCanvas->cd(2);
      gPad->Clear();
      auto *g2 = new TGraph(cx.size(), cz.data(), cx.data());
      g2->SetTitle("XZ;Z_{lab} [mm];X [mm]");
      g2->SetMarkerStyle(20);
      g2->SetMarkerSize(0.6);
      g2->SetMarkerColor(kBlue);
      g2->Draw("AP");
      auto *g2f = new TGraph(sx.size(), sz.data(), sx.data());
      g2f->SetMarkerStyle(24);
      g2f->SetMarkerSize(0.5);
      g2f->SetMarkerColor(kRed);
      g2f->Draw("P SAME");
      drawLegend(g2, g2f);

      fTrackCanvas->cd(3);
      gPad->Clear();
      auto *g3 = new TGraph(cy.size(), cz.data(), cy.data());
      g3->SetTitle("YZ;Z_{lab} [mm];Y [mm]");
      g3->SetMarkerStyle(20);
      g3->SetMarkerSize(0.6);
      g3->SetMarkerColor(kBlue);
      g3->Draw("AP");
      auto *g3f = new TGraph(sy.size(), sz.data(), sy.data());
      g3f->SetMarkerStyle(24);
      g3f->SetMarkerSize(0.5);
      g3f->SetMarkerColor(kRed);
      g3f->Draw("P SAME");
      drawLegend(g3, g3f);

      fTrackCanvas->Update();
   }

   // Update diagnostics
   UpdateDiagnostics(tracks[bestTrack], fittedTracks[0].get());

   // Update info label
   if (fEventInfoLabel)
      fEventInfoLabel->SetText(Form("Ev %d: FITTED p=%.1f MeV/c, KE=%.2f MeV",
                                    fCurrentEvent, pFit, kin.kineticEnergy));
}

// ===========================================================================
// GUI callbacks
// ===========================================================================
void AtUKFDisplay::GuiGotoEvent()
{
   if (fEventEntry)
      GotoEvent(fEventEntry->GetIntNumber());
}

void AtUKFDisplay::GuiNextEvent()
{
   NextEvent();
}

void AtUKFDisplay::GuiPrevEvent()
{
   PrevEvent();
}

void AtUKFDisplay::GuiFit()
{
   FitCurrentTrack();
}

void AtUKFDisplay::GuiRecluster()
{
   if (!fTreeDigi)
      return;

   fTreeDigi->GetEntry(fCurrentEvent);
   if (!fPatEvtArr || fPatEvtArr->GetEntries() == 0)
      return;

   auto *patEvt = (AtPatternEvent *)fPatEvtArr->At(0);
   auto &tracks = patEvt->GetTrackCand();
   if (tracks.empty())
      return;

   // Find largest track
   int bestTrack = 0;
   for (size_t t = 1; t < tracks.size(); t++) {
      if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
         bestTrack = t;
   }

   // Re-cluster with GUI parameters
   int method = fClusterMethodBox ? fClusterMethodBox->GetSelected() : 0;

   tracks[bestTrack].ResetHitClusterArray();
   AtTools::AtTrackTransformer transformer;

   if (method == 0) {
      // Smooth3D
      double radius = fClusterRadiusEntry ? fClusterRadiusEntry->GetNumber() : 10.0;
      double distance = fClusterDistEntry ? fClusterDistEntry->GetNumber() : 20.0;
      transformer.ClusterizeSmooth3D(tracks[bestTrack], radius, distance);
      std::cout << "Smooth3D: radius=" << radius << " distance=" << distance;
   } else {
      // GroupByN
      int hitsPerCluster = fHitsPerClusterEntry ? fHitsPerClusterEntry->GetIntNumber() : 15;
      transformer.ClusterizeByGroup(tracks[bestTrack], hitsPerCluster);
      std::cout << "GroupByN: hitsPerCluster=" << hitsPerCluster;
   }

   int nNewClusters = tracks[bestTrack].GetHitClusterArray()->size();
   std::cout << " → " << nNewClusters << " clusters" << std::endl;

   // Now fit with the new clusters (don't reload from file!)
   FitCurrentTrack(true);
}
