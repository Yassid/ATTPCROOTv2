/// @file cycleEventsGui.C
/// @brief Button-GUI event stepper for the AtGenfitter output (mirrors PUMA/pi_inspect.C).
///        TGMainFrame with [<< Prev] [Next >>] [Event NNN] [Goto] [Save PNG] and a
///        "gated-only" checkbox; an embedded 3-pad canvas (XY / ZX / ZY). 2D canvas
///        only — works over X11 (no Eve/OpenGL).
///
/// Run interactively (NOT -b), needs an X display:
///   root -l 'pipeline/cycleEventsGui.C("/tmp/run_0106_genfitter_gated.root")'
///   root -l 'pipeline/cycleEventsGui.C("/tmp/run_0106_genfitter_gated.root","pid/deuteron_band.json")'
///
/// Gated tracks are drawn solid+colored (vertex star + kinematic stub in red); other
/// tracks are dimmed grey. The header label shows the current fit's theta/KE/chi2/conv.

#include <iostream>
#include <map>
#include <vector>

class FitInspector : public TObject {
public:
   FitInspector(TString file, TString gateFile, Double_t zPadPlane)
   {
      fZpad = zPadPlane;
      fGate = AtTools::AtParticleID::LoadJSON(gateFile.Data());
      if (!fGate.GetCut().IsValid())
         std::cerr << "\n*** WARNING: deuteron gate '" << gateFile << "' did NOT load (empty/invalid).\n"
                   << "***          'gated-only' will show NO events. Launch ROOT from the .../UKF/ dir,\n"
                   << "***          or pass an absolute path as the 2nd argument.\n\n";
      fF = TFile::Open(file);
      if (!fF || fF->IsZombie()) { std::cerr << "cannot open " << file << "\n"; return; }
      fT = (TTree *)fF->Get("cbmsim");
      if (!fT) { std::cerr << "no cbmsim tree\n"; return; }
      fT->SetBranchAddress("AtTrackingEvent", &fTrk);
      fT->SetBranchAddress("AtPIDEvent", &fPid);
      fN = fT->GetEntries();
      RebuildList();
      MakeGui();
      if (!fList.empty())
         Goto(0);
      else {
         if (fLabel) fLabel->SetText(fGate.GetCut().IsValid() ? "no gated events in this file"
                                                              : "gate failed to load -- untick 'gated-only'");
         fCanvas->cd(); fCanvas->Update();
         std::cerr << "step list empty: nothing to draw in gated-only mode (untick the checkbox to see all events)\n";
      }
   }

   void Next() { Goto(fIdx + 1); }
   void Prev() { Goto(fIdx - 1); }
   void Save()
   {
      if (fList.empty()) return;
      TString png = Form("/tmp/cycle_ev%05d.png", fList[fIdx]);
      fCanvas->SaveAs(png);
      std::cout << "saved " << png << "\n";
   }
   void GotoFromEntry()
   {
      if (!fEntry) return;
      int want = (int)fEntry->GetNumberEntry()->GetIntNumber();
      for (size_t k = 0; k < fList.size(); ++k)
         if (fList[k] == want) { Goto((int)k); return; }
      std::cout << "ev " << want << " not in step list (" << (fGatedOnly ? "gated-only" : "all") << ")\n";
   }
   void ToggleGated()
   {
      fGatedOnly = (fCheck && fCheck->IsDown());
      int curEv = fList.empty() ? 0 : fList[fIdx];
      RebuildList();
      // keep position near the event we were on
      fIdx = 0;
      for (size_t k = 0; k < fList.size(); ++k) if (fList[k] >= curEv) { fIdx = (int)k; break; }
      if (!fList.empty()) Goto(fIdx);
      else { fCanvas->Clear(); fCanvas->Update(); if (fLabel) fLabel->SetText("no events in list"); }
   }
   void Goto(int k) // k = index into fList
   {
      if (fList.empty()) return;
      if (k < 0) k = 0;
      if (k >= (int)fList.size()) k = (int)fList.size() - 1;
      fIdx = k;
      if (fEntry) fEntry->GetNumberEntry()->SetIntNumber(fList[fIdx]);
      Draw(fList[fIdx]);
   }

private:
   void RebuildList()
   {
      fList.clear();
      for (Long64_t i = 0; i < fN; ++i) {
         fT->GetEntry(i);
         auto *te = (AtTrackingEvent *)fTrk->At(0);
         if (!te) continue;
         if (!fGatedOnly) { if (!te->GetFittedTracks().empty()) fList.push_back((Int_t)i); continue; }
         auto *pe = (fPid && fPid->GetEntriesFast() > 0) ? (AtPIDEvent *)fPid->At(0) : nullptr;
         bool any = false;
         if (pe) for (auto &r : pe->GetSpyral()) if (r.valid && fGate.IsInside(r.sqrtdEdx, r.brho)) { any = true; break; }
         if (any) fList.push_back((Int_t)i);
      }
      std::cout << "step list: " << fList.size() << " events (" << (fGatedOnly ? "gated-only" : "all-with-fits") << ")\n";
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1850, 720);
      main->SetWindowName("AtGenfitter event inspector");

      auto *bar = new TGHorizontalFrame(main);
      auto *bPrev = new TGTextButton(bar, "<< Prev");
      bPrev->Connect("Clicked()", "FitInspector", this, "Prev()");
      bar->AddFrame(bPrev, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bNext = new TGTextButton(bar, "Next >>");
      bNext->Connect("Clicked()", "FitInspector", this, "Next()");
      bar->AddFrame(bNext, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *lab = new TGLabel(bar, "Event:");
      bar->AddFrame(lab, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 4, 4));
      fEntry = new TGNumberEntry(bar, 0, 7, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0, 1e9);
      bar->AddFrame(fEntry, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bGoto = new TGTextButton(bar, "Goto");
      bGoto->Connect("Clicked()", "FitInspector", this, "GotoFromEntry()");
      bar->AddFrame(bGoto, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *bSave = new TGTextButton(bar, "Save PNG");
      bSave->Connect("Clicked()", "FitInspector", this, "Save()");
      bar->AddFrame(bSave, new TGLayoutHints(kLHintsLeft, 12, 4, 4, 4));

      fCheck = new TGCheckButton(bar, "gated-only");
      fCheck->SetState(fGatedOnly ? kButtonDown : kButtonUp);
      fCheck->Connect("Clicked()", "FitInspector", this, "ToggleGated()");
      bar->AddFrame(fCheck, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 4, 4));

      fLabel = new TGLabel(bar, "                                                            ");
      bar->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 16, 4, 4, 4));

      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *embedded = new TRootEmbeddedCanvas("ec_fitinspect", main, 1830, 660);
      main->AddFrame(embedded, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = embedded->GetCanvas();
      fCanvas->Divide(3, 1); // divide ONCE; Draw() only refreshes the pads

      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   void Draw(int iev)
   {
      fT->GetEntry(iev);
      auto *te = (AtTrackingEvent *)fTrk->At(0);
      auto *pe = (fPid && fPid->GetEntriesFast() > 0) ? (AtPIDEvent *)fPid->At(0) : nullptr;
      if (!te) return;
      auto tracks = te->GetTrackArray();
      auto &fits = te->GetFittedTracks();
      std::map<int, bool> inGate;
      if (pe) for (auto &r : pe->GetSpyral()) inGate[r.trackID] = r.valid && fGate.IsInside(r.sqrtdEdx, r.brho);
      int nGated = 0; for (auto &kv : inGate) if (kv.second) nGated++;

      // header label = representative gated fit's quality
      TString info = Form("ev %d  |  %zu tracks, %zu fits, %d gated", iev, tracks.size(), fits.size(), nGated);
      for (auto &fp : fits) {
         if (!fp) continue;
         if (inGate.count(fp->GetTrackID()) && inGate[fp->GetTrackID()]) {
            const auto &kin = fp->GetKinematics(0);
            auto &m = fp->GetTrackMetadata();
            double c2 = (m && m->GetNdf() > 0) ? m->GetChi2() / m->GetNdf() : -1;
            info += Form("   |   DEUT [%s]: KE=%.1f MeV  #theta=%.1f#circ  #chi^{2}/ndf=%.2f  conv=%d",
                         (m && m->GetGoodFit()) ? "GOOD" : "FLAGGED", kin.kineticEnergy, kin.theta * TMath::RadToDeg(),
                         c2, m ? m->GetFitConverged() : -1);
            break;
         }
      }
      if (fLabel) fLabel->SetText(info);

      const Color_t pal[] = {kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1, kCyan + 2, kMagenta + 1, kAzure + 7};
      auto pick = [&](int proj, double X, double Y, double Z, double &a, double &b) {
         if (proj == 0) { a = X; b = Y; } else if (proj == 1) { a = Z; b = X; } else { a = Z; b = Y; }
      };
      const char *axX[3] = {"X [mm]", "z_{lab} [mm]", "z_{lab} [mm]"};
      const char *axY[3] = {"Y [mm]", "X [mm]", "Y [mm]"};

      for (int proj = 0; proj < 3; ++proj) {
         fCanvas->cd(proj + 1); gPad->Clear(); gPad->SetGrid();
         double aLo, aHi, bLo, bHi;
         if (proj == 0) { aLo = bLo = -250; aHi = bHi = 250; }
         else { aLo = -50; aHi = fZpad + 50; bLo = -250; bHi = 250; }
         gPad->DrawFrame(aLo, bLo, aHi, bHi, Form("ev %d  %s vs %s;%s;%s", iev, axY[proj], axX[proj], axX[proj], axY[proj]));
         if (proj == 0) {
            auto *lg = new TLatex(); lg->SetNDC(); lg->SetTextSize(0.033);
            lg->SetTextColor(kBlue + 1); lg->DrawLatex(0.15, 0.86, "#bullet hits");
            lg->SetTextColor(kBlack);    lg->DrawLatex(0.15, 0.82, "#Box clusters (fit input)");
            lg->SetTextColor(kGreen + 2); lg->DrawLatex(0.15, 0.78, "#minus fit GOOD");
            lg->SetTextColor(kOrange + 7); lg->DrawLatex(0.15, 0.74, "#minus fit FLAGGED   #star vertex");
         }
         for (size_t i = 0; i < tracks.size(); ++i) {
            auto &tr = tracks[i];
            bool g = inGate.count(tr.GetTrackID()) ? inGate[tr.GetTrackID()] : false;
            auto &hits = tr.GetHitArray();
            if (hits.empty()) continue;
            auto *gr = new TGraph();
            for (auto &h : hits) { const auto &p = h->GetPosition(); double a, b; pick(proj, p.X(), p.Y(), fZpad - p.Z(), a, b); gr->SetPoint(gr->GetN(), a, b); }
            gr->SetMarkerStyle(g ? 20 : 24); gr->SetMarkerSize(g ? 0.55 : 0.4);
            gr->SetMarkerColor(g ? pal[i % 7] : (kGray + 1)); gr->Draw("P SAME");

            // clusters the fit actually consumed (black open squares), gated tracks only
            if (g) {
               auto *hcArr = tr.GetHitClusterArray();
               if (hcArr && !hcArr->empty()) {
                  auto *gc = new TGraph();
                  for (auto &cl : *hcArr) { auto cp = cl.GetPosition(); double a, b; pick(proj, cp.X(), cp.Y(), fZpad - cp.Z(), a, b); gc->SetPoint(gc->GetN(), a, b); }
                  gc->SetMarkerStyle(25); gc->SetMarkerSize(0.9); gc->SetMarkerColor(kBlack); gc->Draw("P SAME");
               }
            }
         }
         for (auto &fp : fits) {
            if (!fp) continue;
            bool g = inGate.count(fp->GetTrackID()) ? inGate[fp->GetTrackID()] : false;
            const auto &kin = fp->GetKinematics(0); const auto &v = fp->GetVertex(0);
            auto &mq = fp->GetTrackMetadata();
            bool good = mq && mq->GetGoodFit();
            Color_t cc = !g ? (kGray + 2) : (good ? kGreen + 2 : kOrange + 7); // green=good, orange=flagged
            double a, b; pick(proj, v.X(), v.Y(), v.Z(), a, b);
            auto *m = new TMarker(a, b, 29); m->SetMarkerSize(g ? 2.0 : 1.2); m->SetMarkerColor(cc); m->Draw();

            // the fitted trajectory (genfit smoothed points, already lab-frame mm -> no z transform)
            const auto &traj = fp->GetSmoothedPositions();
            if (!traj.empty()) {
               auto *gt = new TGraph();
               for (auto &sp : traj) { double a1, b1; pick(proj, sp.X(), sp.Y(), sp.Z(), a1, b1); gt->SetPoint(gt->GetN(), a1, b1); }
               gt->SetLineColor(cc); gt->SetLineWidth(2); gt->Draw("L SAME");
            } else if (kin.kineticEnergy > 0 && std::isfinite(kin.theta) && std::isfinite(kin.phi)) {
               double L = 60, ux = std::sin(kin.theta) * std::cos(kin.phi), uy = std::sin(kin.theta) * std::sin(kin.phi), uz = std::cos(kin.theta);
               double a1, b1; pick(proj, v.X() + L * ux, v.Y() + L * uy, v.Z() + L * uz, a1, b1);
               auto *l = new TLine(a, b, a1, b1); l->SetLineColor(cc); l->SetLineWidth(g ? 3 : 1); l->SetLineStyle(2); l->Draw();
            }
         }
      }
      fCanvas->Modified(); fCanvas->Update();
      gSystem->ProcessEvents();
   }

   TFile *fF{nullptr};
   TTree *fT{nullptr};
   TClonesArray *fTrk{nullptr}, *fPid{nullptr};
   AtTools::AtParticleID fGate;
   Double_t fZpad{1000.0};
   Long64_t fN{0};
   std::vector<Int_t> fList;
   int fIdx{0};
   Bool_t fGatedOnly{kTRUE};
   TCanvas *fCanvas{nullptr};
   TGNumberEntry *fEntry{nullptr};
   TGLabel *fLabel{nullptr};
   TGCheckButton *fCheck{nullptr};

   ClassDef(FitInspector, 0);
};

void cycleEventsGui(TString file = "/tmp/run_0106_genfitter_gated.root",
                    TString gateFile =
                       "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF/pid/deuteron_band.json",
                    Double_t zPadPlane = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   new FitInspector(file, gateFile, zPadPlane);
}
