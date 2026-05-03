/// @file pi_inspect.C
/// @brief Interactive event-by-event viewer for the π+/π- fit projections.
///
/// Wraps inspect_fit.C's drawing logic in a small TGMainFrame with
/// Prev / Next / Goto controls and an embedded canvas. Per event:
///   gray dots     — all MC truth points (cm → mm)
///   coloured dots — PRA hits, one colour per track
///   red ⋆        — fitted vertex (initialPositionXtr)
///   red triangle  — first cluster (initialPosition)
///   red ◦ + line  — UKF smoothed positions (line through them)
///   red dashed    — kinematic-direction stub from vertex (50 mm)
///   black ⋆      — MC truth primary vertex
///   black dashed  — MC primary direction (50 mm)
///
/// Three pads (XY, XZ, YZ) plus a fourth info pad with fit/truth tables.
///
/// Run: root -l pi_inspect.C
///      root -l 'pi_inspect.C("./data/output_digi.root", "./data/attpcsim.root")'

#include <iostream>
#include <vector>

// TObject base + ClassDef are needed for ROOT's signal/slot Connect to find
// the class's method table when the macro is interpreted by Cling.
class PiInspector : public TObject {
public:
   PiInspector(TString digiFile, TString simFile)
   {
      fD = TFile::Open(digiFile);
      fS = TFile::Open(simFile);
      if (!fD || fD->IsZombie() || !fS || fS->IsZombie()) {
         std::cerr << "Cannot open inputs\n";
         return;
      }
      fT_D = (TTree *)fD->Get("cbmsim");
      fT_S = (TTree *)fS->Get("cbmsim");
      fMcPts = new TClonesArray("AtMCPoint");
      fMcTrks = new TClonesArray("AtMCTrack");
      fT_D->SetBranchAddress("AtPatternEvent", &fPeArr);
      if (fT_D->GetBranch("AtTrackingEvent"))
         fT_D->SetBranchAddress("AtTrackingEvent", &fTeArr);
      fT_S->SetBranchAddress("AtTpcPoint", &fMcPts);
      fT_S->SetBranchAddress("MCTrack", &fMcTrks);
      fNEvents = std::min(fT_D->GetEntries(), fT_S->GetEntries());
      std::cout << "PiInspector: " << fNEvents << " events available\n";
      MakeGui();
      Draw(0);
   }

   void Next() { Goto(fCurrent + 1); }
   void Prev() { Goto(fCurrent - 1); }
   void GotoFromEntry()
   {
      if (fEntry)
         Goto((int)fEntry->GetNumberEntry()->GetIntNumber());
   }
   void Goto(int e)
   {
      if (e < 0)
         e = 0;
      if (e >= fNEvents)
         e = fNEvents - 1;
      fCurrent = e;
      if (fEntry)
         fEntry->GetNumberEntry()->SetIntNumber(fCurrent);
      Draw(fCurrent);
   }

private:
   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1300, 750);
      main->SetWindowName("PiInspector");

      // Top control bar: [Prev] [Next] [Event NNN] [Goto]
      auto *bar = new TGHorizontalFrame(main);
      auto *bPrev = new TGTextButton(bar, "<< Prev");
      bPrev->Connect("Clicked()", "PiInspector", this, "Prev()");
      bar->AddFrame(bPrev, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *bNext = new TGTextButton(bar, "Next >>");
      bNext->Connect("Clicked()", "PiInspector", this, "Next()");
      bar->AddFrame(bNext, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *lab = new TGLabel(bar, "Event:");
      bar->AddFrame(lab, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 4, 4));

      fEntry = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0, 1e9);
      bar->AddFrame(fEntry, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *bGoto = new TGTextButton(bar, "Goto");
      bGoto->Connect("Clicked()", "PiInspector", this, "GotoFromEntry()");
      bar->AddFrame(bGoto, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      // Embedded canvas (split into 2x2)
      auto *embedded = new TRootEmbeddedCanvas("ec_pi_inspect", main, 1280, 700);
      main->AddFrame(embedded, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = embedded->GetCanvas();
      fCanvas->Divide(2, 2);

      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   void Draw(int event)
   {
      if (event < 0 || event >= fNEvents)
         return;
      fT_D->GetEntry(event);
      fT_S->GetEntry(event);

      auto *pe = (AtPatternEvent *)(fPeArr ? fPeArr->At(0) : nullptr);
      const std::vector<AtTrack> *tracks = pe ? &pe->GetTrackCand() : nullptr;

      AtTrackingEvent *te = nullptr;
      if (fTeArr && fTeArr->GetEntries() > 0)
         te = (AtTrackingEvent *)fTeArr->At(0);

      // --- Three projection pads ---
      DrawProjection(1, 0, event, tracks, te);
      DrawProjection(2, 1, event, tracks, te);
      DrawProjection(3, 2, event, tracks, te);

      // --- Info pad ---
      fCanvas->cd(4);
      gPad->Clear();
      auto *pt = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
      pt->SetTextAlign(12);
      pt->SetTextFont(102);
      pt->SetTextSize(0.045);
      pt->SetFillColor(0);
      pt->SetBorderSize(0);
      pt->AddText(Form("Event %d / %d", event, (int)fNEvents - 1));
      pt->AddText(Form("MC primaries: %d   PRA tracks: %d   UKF fits: %d",
                       fMcTrks->GetEntries(),
                       tracks ? (int)tracks->size() : 0,
                       te ? (int)te->GetFittedTracks().size() : 0));
      pt->AddText(" ");
      pt->AddLine(0., 0., 1., 0.);

      pt->AddText("--- Truth primaries ---");
      for (int k = 0; k < fMcTrks->GetEntries() && k < 6; ++k) {
         auto *m = (AtMCTrack *)fMcTrks->At(k);
         int pdg = m->GetPdgCode();
         double pmag = std::sqrt(m->GetPx() * m->GetPx() + m->GetPy() * m->GetPy()
                                 + m->GetPz() * m->GetPz());
         if (pmag <= 0)
            continue;
         double th = std::acos(m->GetPz() / pmag) * 180.0 / TMath::Pi();
         double ph = std::atan2(m->GetPy(), m->GetPx()) * 180.0 / TMath::Pi();
         pt->AddText(Form("  pdg=%d  |p|=%.0f MeV  theta=%.1f  phi=%.1f", pdg, pmag * 1000, th, ph));
      }

      // Cache MC point (x,y) and trackID for the truth-match below.
      const int nMC = fMcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC);
      std::vector<int> mcTid(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)fMcPts->At(k);
         mcX[k] = mp->GetX() * 10.0;
         mcY[k] = mp->GetY() * 10.0;
         mcTid[k] = mp->GetTrackID();
      }
      auto truthMatch = [&](const AtTrack &tr) -> std::string {
         std::map<int, int> votes;
         const double tol2 = 9.0; // 3 mm
         for (auto &hp : tr.GetHitArray()) {
            const auto &p = hp->GetPosition();
            double bestD2 = tol2;
            int bestTid = -1;
            for (int k = 0; k < nMC; ++k) {
               double dx = p.X() - mcX[k];
               double dy = p.Y() - mcY[k];
               double d2 = dx * dx + dy * dy;
               if (d2 < bestD2) { bestD2 = d2; bestTid = mcTid[k]; }
            }
            if (bestTid >= 0) votes[bestTid]++;
         }
         if (votes.empty()) return "no-MC";
         int bestTid = -1, bestN = 0;
         for (auto &kv : votes) if (kv.second > bestN) { bestN = kv.second; bestTid = kv.first; }
         if (bestTid >= 0 && bestTid < fMcTrks->GetEntries()) {
            int pdg = ((AtMCTrack *)fMcTrks->At(bestTid))->GetPdgCode();
            return Form("tid=%d pdg=%d", bestTid, pdg);
         }
         return Form("tid=%d (sec)", bestTid);
      };

      pt->AddText(" ");
      pt->AddText(Form("--- PRA tracks (%d) ---", tracks ? (int)tracks->size() : 0));
      if (tracks) {
         for (size_t i = 0; i < tracks->size() && i < 8; ++i) {
            std::string tag = truthMatch((*tracks)[i]);
            pt->AddText(Form("  Tr %zu  hits=%-4d  %s", i,
                             (int)(*tracks)[i].GetHitArray().size(), tag.c_str()));
         }
      }

      pt->AddText(" ");
      pt->AddText("--- UKF fits ---");
      if (te) {
         const auto &fits = te->GetFittedTracks();
         for (size_t i = 0; i < fits.size() && i < 6; ++i) {
            const auto &ft = *fits[i];
            const auto &kin = ft.GetKinematics();
            const auto &pi = ft.GetParticleInfo();
            pt->AddText(Form("  pdg=%-4s  KE=%.1f MeV  theta=%.1f  phi=%.1f",
                             pi.idPDG.Data(), kin.kineticEnergy,
                             kin.theta * 180.0 / TMath::Pi(),
                             kin.phi * 180.0 / TMath::Pi()));
         }
      }
      pt->Draw();
      gPad->Modified();
      gPad->Update();

      fCanvas->Update();
   }

   void DrawProjection(int padNum, int proj, int event,
                       const std::vector<AtTrack> *tracks, AtTrackingEvent *te)
   {
      fCanvas->cd(padNum);
      gPad->Clear();
      gPad->SetGrid();

      double aLo, aHi, bLo, bHi;
      const char *axA, *axB;
      if (proj == 0) {
         aLo = -150; aHi = 150; bLo = -150; bHi = 150;
         axA = "X [mm]"; axB = "Y [mm]";
      } else if (proj == 1) {
         aLo = -50; aHi = 350; bLo = -150; bHi = 150;
         axA = "Z [mm]"; axB = "X [mm]";
      } else {
         aLo = -50; aHi = 350; bLo = -150; bHi = 150;
         axA = "Z [mm]"; axB = "Y [mm]";
      }

      gPad->DrawFrame(aLo, bLo, aHi, bHi,
                      Form("event %d  %s vs %s ;%s;%s", event, axB, axA, axA, axB));

      auto pick = [&](double X, double Y, double Z, double &a, double &b) {
         if (proj == 0) { a = X; b = Y; }
         else if (proj == 1) { a = Z; b = X; }
         else                { a = Z; b = Y; }
      };

      // MC truth points, colored per PDG of their parent MCTrack.
      // π+ → green, π− → magenta, anything else (secondaries / non-pion
      // primaries) → gray.
      auto *gPiPlus = new TGraph();
      auto *gPiMinus = new TGraph();
      auto *gOther = new TGraph();
      for (int k = 0; k < fMcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)fMcPts->At(k);
         int tid = mp->GetTrackID();
         int pdg = 0;
         if (tid >= 0 && tid < fMcTrks->GetEntries())
            pdg = ((AtMCTrack *)fMcTrks->At(tid))->GetPdgCode();
         double aP, bP;
         pick(mp->GetX() * 10.0, mp->GetY() * 10.0, mp->GetZ() * 10.0, aP, bP);
         if (pdg == 211) gPiPlus->SetPoint(gPiPlus->GetN(), aP, bP);
         else if (pdg == -211) gPiMinus->SetPoint(gPiMinus->GetN(), aP, bP);
         else gOther->SetPoint(gOther->GetN(), aP, bP);
      }
      gOther->SetMarkerStyle(7);
      gOther->SetMarkerColor(kGray + 1);
      gOther->Draw("P SAME");
      gPiPlus->SetMarkerStyle(7);
      gPiPlus->SetMarkerColor(kGreen + 2);
      gPiPlus->Draw("P SAME");
      gPiMinus->SetMarkerStyle(7);
      gPiMinus->SetMarkerColor(kMagenta + 1);
      gPiMinus->Draw("P SAME");

      // PRA tracks
      const Color_t cols[] = {kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1,
                              kCyan + 2, kYellow - 2, kAzure + 7};
      if (tracks) {
         for (size_t i = 0; i < tracks->size(); ++i) {
            auto &h = (*tracks)[i].GetHitArray();
            if (h.empty())
               continue;
            auto *g = new TGraph();
            for (auto &hp : h) {
               const auto &p = hp->GetPosition();
               double a, b;
               pick(p.X(), p.Y(), p.Z(), a, b);
               g->SetPoint(g->GetN(), a, b);
            }
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.5);
            g->SetMarkerColor(cols[i % 7]);
            g->Draw("P SAME");
         }
      }

      // Fitted tracks
      if (te) {
         const auto &fits = te->GetFittedTracks();
         for (size_t i = 0; i < fits.size(); ++i) {
            const auto &ft = *fits[i];
            const auto &props = ft.GetTrackPropertiesStruct();
            const auto &kin = ft.GetKinematics();
            const auto &pi = ft.GetParticleInfo();
            const auto &pV = props.initialPositionXtr;
            const auto &p0 = props.initialPosition;

            // vertex (red ⋆)
            if (pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0) {
               double a, b;
               pick(pV.X(), pV.Y(), pV.Z(), a, b);
               auto *m = new TMarker(a, b, 29);
               m->SetMarkerSize(2.0);
               m->SetMarkerColor(kRed + 1);
               m->Draw();
            }
            // first cluster (red ▲)
            if (p0.X() != 0 || p0.Y() != 0 || p0.Z() != 0) {
               double a, b;
               pick(p0.X(), p0.Y(), p0.Z(), a, b);
               auto *m = new TMarker(a, b, 22);
               m->SetMarkerSize(1.2);
               m->SetMarkerColor(kRed + 1);
               m->Draw();
            }
            // smoothed UKF positions (red ◦ + line)
            const auto &smoothed = props.fSmoothedPositions;
            if (!smoothed.empty()) {
               auto *gSm = new TGraph();
               for (const auto &sp : smoothed) {
                  double a, b;
                  pick(sp.X(), sp.Y(), sp.Z(), a, b);
                  gSm->SetPoint(gSm->GetN(), a, b);
               }
               gSm->SetMarkerStyle(24);
               gSm->SetMarkerSize(0.7);
               gSm->SetLineColor(kRed + 1);
               gSm->SetLineWidth(2);
               gSm->SetMarkerColor(kRed + 1);
               gSm->Draw("LP SAME");
            }
            // kinematic-direction stub (red dashed)
            if (kin.kineticEnergy > 0 && pi.mass > 0
                && std::isfinite(kin.theta) && std::isfinite(kin.phi)
                && (pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0)) {
               const double L = 50.0;
               const double ux = std::sin(kin.theta) * std::cos(kin.phi);
               const double uy = std::sin(kin.theta) * std::sin(kin.phi);
               const double uz = std::cos(kin.theta);
               double a0, b0, a1, b1;
               pick(pV.X(), pV.Y(), pV.Z(), a0, b0);
               pick(pV.X() + L * ux, pV.Y() + L * uy, pV.Z() + L * uz, a1, b1);
               auto *l = new TLine(a0, b0, a1, b1);
               l->SetLineColor(kRed + 1);
               l->SetLineWidth(2);
               l->SetLineStyle(2);
               l->Draw();
            }
         }
      }

      // MC primaries (black ⋆ + dashed direction)
      for (int k = 0; k < std::min(6, fMcTrks->GetEntries()); ++k) {
         auto *mt = (AtMCTrack *)fMcTrks->At(k);
         int pdg = mt->GetPdgCode();
         if (std::abs(pdg) != 211 && std::abs(pdg) != 321 && std::abs(pdg) != 2212)
            continue;
         double vx = mt->GetStartX() * 10.0;
         double vy = mt->GetStartY() * 10.0;
         double vz = mt->GetStartZ() * 10.0;
         double a, b;
         pick(vx, vy, vz, a, b);
         auto *m = new TMarker(a, b, 29);
         m->SetMarkerSize(2.0);
         m->SetMarkerColor(kBlack);
         m->Draw();
         double px = mt->GetPx(), py = mt->GetPy(), pz = mt->GetPz();
         double pmag = std::sqrt(px * px + py * py + pz * pz);
         if (pmag > 0) {
            const double L = 50.0;
            double a1, b1;
            pick(vx + L * px / pmag, vy + L * py / pmag, vz + L * pz / pmag, a1, b1);
            auto *l = new TLine(a, b, a1, b1);
            l->SetLineColor(kBlack);
            l->SetLineWidth(2);
            l->SetLineStyle(3);
            l->Draw();
         }
      }

      // Compact legend (top-right). Dummy graphs/markers/lines just for
      // the legend symbols — they're not drawn anywhere themselves because
      // their TGraph has no points and the markers/lines are not Drawn().
      auto *leg = new TLegend(0.65, 0.55, 0.985, 0.92);
      leg->SetTextSize(0.028);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      auto *lgPiP = new TGraph(); lgPiP->SetMarkerStyle(7); lgPiP->SetMarkerColor(kGreen + 2);
      auto *lgPiM = new TGraph(); lgPiM->SetMarkerStyle(7); lgPiM->SetMarkerColor(kMagenta + 1);
      auto *lgSec = new TGraph(); lgSec->SetMarkerStyle(7); lgSec->SetMarkerColor(kGray + 1);
      auto *lgPRA = new TGraph(); lgPRA->SetMarkerStyle(20); lgPRA->SetMarkerSize(0.6); lgPRA->SetMarkerColor(kBlue + 1);
      auto *lgSm  = new TGraph(); lgSm->SetMarkerStyle(24); lgSm->SetMarkerSize(0.7); lgSm->SetLineColor(kRed + 1); lgSm->SetMarkerColor(kRed + 1);
      auto *lgVtxFit = new TMarker(0, 0, 29); lgVtxFit->SetMarkerColor(kRed + 1);
      auto *lgFirst  = new TMarker(0, 0, 22); lgFirst->SetMarkerColor(kRed + 1);
      auto *lgVtxMC  = new TMarker(0, 0, 29); lgVtxMC->SetMarkerColor(kBlack);
      auto *lgDirFit = new TLine(); lgDirFit->SetLineColor(kRed + 1); lgDirFit->SetLineStyle(2); lgDirFit->SetLineWidth(2);
      auto *lgDirMC  = new TLine(); lgDirMC->SetLineColor(kBlack); lgDirMC->SetLineStyle(3); lgDirMC->SetLineWidth(2);
      leg->AddEntry(lgPiP,    "MC #pi+ point",        "P");
      leg->AddEntry(lgPiM,    "MC #pi- point",        "P");
      leg->AddEntry(lgSec,    "MC secondary",          "P");
      leg->AddEntry(lgPRA,    "PRA hit (per track)",   "P");
      leg->AddEntry(lgSm,     "UKF smoothed",          "LP");
      leg->AddEntry(lgVtxFit, "UKF vertex (Xtr)",      "P");
      leg->AddEntry(lgFirst,  "UKF first cluster",     "P");
      leg->AddEntry(lgDirFit, "UKF kin. direction",    "L");
      leg->AddEntry(lgVtxMC,  "MC primary vertex",     "P");
      leg->AddEntry(lgDirMC,  "MC primary direction",  "L");
      leg->Draw();

      gPad->Modified();
      gPad->Update();
   }

   TFile *fD{nullptr}, *fS{nullptr};
   TTree *fT_D{nullptr}, *fT_S{nullptr};
   TClonesArray *fPeArr{nullptr}, *fTeArr{nullptr};
   TClonesArray *fMcPts{nullptr}, *fMcTrks{nullptr};
   Long64_t fNEvents{0};
   int fCurrent{0};
   TCanvas *fCanvas{nullptr};
   TGNumberEntry *fEntry{nullptr};

   ClassDef(PiInspector, 0);
};

void pi_inspect(TString digiFile = "./data/output_digi.root",
                TString simFile = "./data/attpcsim.root")
{
   new PiInspector(digiFile, simFile);
}
