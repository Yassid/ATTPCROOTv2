#include "AtUKFDisplay.h"

#include "AtELossCATIMA.h"
#include "AtFittedTrack.h"
#include "AtFitterUKF.h"
#include "AtHitCluster.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtTrackingEvent.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveStraightLineSet.h>
#include <TFile.h>
#include <TGButton.h>
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
#include <TRootBrowser.h>
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
   if (!gEve)
      TEveManager::Create();

   MakeGui();
   MakeDiagnosticsCanvas();
}

AtUKFDisplay::~AtUKFDisplay()
{
   fInstance = nullptr;
}

// ===========================================================================
// File loading
// ===========================================================================
void AtUKFDisplay::LoadFiles(const char *digiFile, const char *fittedFile)
{
   fFileDigi.reset(TFile::Open(digiFile));
   if (!fFileDigi || fFileDigi->IsZombie()) {
      std::cerr << "AtUKFDisplay: cannot open " << digiFile << std::endl;
      return;
   }
   fTreeDigi = (TTree *)fFileDigi->Get("cbmsim");
   fTreeDigi->SetBranchAddress("AtPatternEvent", &fPatEvtArr);
   fNEvents = fTreeDigi->GetEntries();

   if (fittedFile) {
      fFileFit.reset(TFile::Open(fittedFile));
      if (fFileFit && !fFileFit->IsZombie()) {
         fTreeFit = (TTree *)fFileFit->Get("cbmsim");
         fTreeFit->SetBranchAddress("AtTrackingEvent", &fTrackingEvtArr);
      }
   }

   std::cout << "AtUKFDisplay: loaded " << fNEvents << " events from " << digiFile << std::endl;
   if (fTreeFit)
      std::cout << "AtUKFDisplay: fitted data loaded from " << fittedFile << std::endl;
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
   auto *browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   fMainFrame = new TGMainFrame(gClient->GetRoot(), 250, 600);
   fMainFrame->SetWindowName("UKF Fitter");
   fMainFrame->SetCleanup(kDeepCleanup);

   MakeControlPanel(fMainFrame);

   fMainFrame->MapSubwindows();
   fMainFrame->Resize();
   fMainFrame->MapWindow();

   browser->StopEmbedding("UKF Fitter");
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

   // Checkboxes
   fStragglingBtn = new TGCheckButton(vf, "Energy straggling");
   vf->AddFrame(fStragglingBtn, new TGLayoutHints(kLHintsLeft, 5, 2, 3, 1));

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
   // Create a diagnostics canvas in the Eve browser's bottom area
   auto *browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kBottom);

   fDiagCanvas = new TCanvas("cUKFDiag", "UKF Diagnostics", 1000, 300);
   fDiagCanvas->Divide(4, 1);

   browser->StopEmbedding("Diagnostics");
}

// ===========================================================================
// Drawing
// ===========================================================================
void AtUKFDisplay::ClearEveElements()
{
   if (fClusterPoints) {
      fClusterPoints->Destroy();
      fClusterPoints = nullptr;
   }
   if (fSmoothedPoints) {
      fSmoothedPoints->Destroy();
      fSmoothedPoints = nullptr;
   }
   if (fFilteredPoints) {
      fFilteredPoints->Destroy();
      fFilteredPoints = nullptr;
   }
   if (fResidualLines) {
      fResidualLines->Destroy();
      fResidualLines = nullptr;
   }
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
      gEve->Redraw3D();
      return;
   }

   auto *patEvt = (AtPatternEvent *)fPatEvtArr->At(0);
   auto &tracks = patEvt->GetTrackCand();
   if (tracks.empty()) {
      if (fEventInfoLabel)
         fEventInfoLabel->SetText(Form("Event %d: no tracks", fCurrentEvent));
      gEve->Redraw3D();
      return;
   }

   // Use the largest track
   int bestTrack = 0;
   for (size_t t = 1; t < tracks.size(); t++) {
      if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
         bestTrack = t;
   }
   auto &track = tracks[bestTrack];

   // Draw clusters
   if (fDrawOptions.find('C') != std::string::npos)
      DrawClusters(track);

   // Draw fitted track if available
   AtFittedTrack *fitted = nullptr;
   if (fTrackingEvtArr && fTrackingEvtArr->GetEntries() > 0) {
      auto *te = (AtTrackingEvent *)fTrackingEvtArr->At(0);
      auto &ft = te->GetFittedTracks();
      if (!ft.empty())
         fitted = ft[0].get();
   }

   if (fitted && fDrawOptions.find('S') != std::string::npos)
      DrawFittedTrack(*fitted);

   // Update info label
   int nClusters = track.GetHitClusterArray()->size();
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

   gEve->Redraw3D(kTRUE);

   // Auto-fit if enabled
   if (fAutoFitBtn && fAutoFitBtn->IsOn() && !fitted)
      FitCurrentTrack();
}

void AtUKFDisplay::DrawClusters(AtTrack &track)
{
   auto *clusters = track.GetHitClusterArray();
   fClusterPoints = new TEvePointSet("Clusters", clusters->size(), TEvePointSelectorConsumer::kTVT_XYZ);
   fClusterPoints->SetMarkerColor(kBlue);
   fClusterPoints->SetMarkerSize(1.5);
   fClusterPoints->SetMarkerStyle(20);

   for (size_t i = 0; i < clusters->size(); i++) {
      auto pos = clusters->at(i).GetPosition();
      // Convert to lab frame for 3D display
      fClusterPoints->SetPoint(i, pos.X(), pos.Y(), fZPadPlane - pos.Z());
   }

   gEve->AddElement(fClusterPoints);
}

void AtUKFDisplay::DrawFittedTrack(const AtFittedTrack &fitted)
{
   auto &smoothed = fitted.GetSmoothedPositions();
   auto vtx = fitted.GetVertex();

   fSmoothedPoints = new TEvePointSet("Smoothed", smoothed.size() + 1, TEvePointSelectorConsumer::kTVT_XYZ);
   fSmoothedPoints->SetMarkerColor(kRed);
   fSmoothedPoints->SetMarkerSize(1.2);
   fSmoothedPoints->SetMarkerStyle(24);

   // Vertex point
   fSmoothedPoints->SetPoint(0, vtx.X(), vtx.Y(), vtx.Z());

   for (size_t i = 0; i < smoothed.size(); i++) {
      fSmoothedPoints->SetPoint(i + 1, smoothed[i].X(), smoothed[i].Y(), smoothed[i].Z());
   }

   gEve->AddElement(fSmoothedPoints);
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
      }
   }

   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1); // TODO: adjust for particle type
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1});
   eloss->SetMaterial(mat);

   fFitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass, std::move(eloss));
   fFitter->SetBField({0, 0, 2.85});

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
   fFitter->SetMomentumSigmaFrac(0.3);

   fZPadPlane = zPadPlane;
}

void AtUKFDisplay::FitCurrentTrack()
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

   // Recreate fitter with current GUI parameters
   CreateFitter();
   fFitter->Init();

   // Fit using the public FitEvent interface
   AtTrackingEvent trackingEvent;
   fFitter->FitEvent(&trackingEvent, patEvt);

   auto &fittedTracks = trackingEvent.GetFittedTracks();
   if (fittedTracks.empty()) {
      std::cout << "AtUKFDisplay: fit failed for event " << fCurrentEvent << std::endl;
      return;
   }

   std::cout << "AtUKFDisplay: fitted event " << fCurrentEvent << " — KE="
             << fittedTracks[0]->GetKinematics().kineticEnergy << " MeV" << std::endl;

   // Redraw with the new fit result
   ClearEveElements();

   // Use largest track
   int bestTrack = 0;
   for (size_t t = 1; t < tracks.size(); t++) {
      if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
         bestTrack = t;
   }

   if (fDrawOptions.find('C') != std::string::npos)
      DrawClusters(tracks[bestTrack]);
   if (fDrawOptions.find('S') != std::string::npos)
      DrawFittedTrack(*fittedTracks[0]);

   UpdateDiagnostics(tracks[bestTrack], fittedTracks[0].get());

   gEve->Redraw3D(kTRUE);
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
