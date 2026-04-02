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
void AtUKFDisplay::LoadFiles(const char *digiFile, const char *fittedFile)
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
      // TCanvas mode: 2D projections
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

      fTrackCanvas->cd(1);
      gPad->Clear();
      auto *gXY = new TGraph(cx.size(), cx.data(), cy.data());
      gXY->SetTitle(Form("Event %d  XY;X [mm];Y [mm]", fCurrentEvent));
      gXY->SetMarkerStyle(20);
      gXY->SetMarkerSize(0.6);
      gXY->SetMarkerColor(kBlue);
      gXY->Draw("AP");
      if (!sx.empty()) {
         auto *gXYf = new TGraph(sx.size(), sx.data(), sy.data());
         gXYf->SetMarkerStyle(24);
         gXYf->SetMarkerSize(0.5);
         gXYf->SetMarkerColor(kRed);
         gXYf->Draw("P SAME");
      }

      fTrackCanvas->cd(2);
      gPad->Clear();
      auto *gXZ = new TGraph(cx.size(), cz.data(), cx.data());
      gXZ->SetTitle("XZ;Z_{lab} [mm];X [mm]");
      gXZ->SetMarkerStyle(20);
      gXZ->SetMarkerSize(0.6);
      gXZ->SetMarkerColor(kBlue);
      gXZ->Draw("AP");
      if (!sz.empty()) {
         auto *gXZf = new TGraph(sx.size(), sz.data(), sx.data());
         gXZf->SetMarkerStyle(24);
         gXZf->SetMarkerSize(0.5);
         gXZf->SetMarkerColor(kRed);
         gXZf->Draw("P SAME");
      }

      fTrackCanvas->cd(3);
      gPad->Clear();
      auto *gYZ = new TGraph(cy.size(), cz.data(), cy.data());
      gYZ->SetTitle("YZ;Z_{lab} [mm];Y [mm]");
      gYZ->SetMarkerStyle(20);
      gYZ->SetMarkerSize(0.6);
      gYZ->SetMarkerColor(kBlue);
      gYZ->Draw("AP");
      if (!sz.empty()) {
         auto *gYZf = new TGraph(sy.size(), sz.data(), sy.data());
         gYZf->SetMarkerStyle(24);
         gYZf->SetMarkerSize(0.5);
         gYZf->SetMarkerColor(kRed);
         gYZf->Draw("P SAME");
      }

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

   // No Eve redraw needed in TCanvas mode
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
