/// @file pra_explorer.C
/// @brief Interactive PRA parameter explorer. Loads a DIGITIZED file (AtEventH = PSA
///        hits), and for the current event re-runs the TriplClust track finder LIVE with
///        whatever parameters you type, drawing the found tracks + their fitted circles
///        on the pad plane. Cycle events with Prev/Next/Goto; edit any parameter and hit
///        "Run PRA" to see the effect immediately. No re-digitization, no Eve (2D X11
///        canvas, WSLg-safe).
///
/// Params exposed (TriplClust: Dalitz et al.): t=Tcluster s=Scluster k=Ktriplet
/// n=Ntriplet m=Mcluster r=Rsmooth a=Atriplet; clustering radius/distance; arc-walk
/// target clusters; circle-merge (annular fragment merge) on/off + tol; SelectAndMerge
/// on/off; chargeFromCenter on/off; min hits/track.
///
/// Run: root -l 'pra_explorer.C("/mnt/f/puma_sweep/output_digi_pi200_primary.root", 200)'
///      (optional 3rd arg: sim file to overlay MC truth points)

#include <iostream>
#include <vector>

// NOTE: do NOT #include the ATTPCROOT headers here — Cling hangs parsing their heavy
// transitive includes. The classes autoload via rootmap once libAtReconstruction is
// loaded in the entry function below (same pattern as pi_inspect.C).

class PRAExplorer : public TObject {
public:
   PRAExplorer(TString digiFile, double p0 = 0.0, TString simFile = "")
   {
      fP0 = p0;
      fD = TFile::Open(digiFile);
      if (!fD || fD->IsZombie()) { std::cerr << "Cannot open " << digiFile << "\n"; return; }
      fT_D = (TTree *)fD->Get("cbmsim");
      fEvArr = new TClonesArray("AtEvent");
      fT_D->SetBranchAddress("AtEventH", &fEvArr);
      fNEvents = fT_D->GetEntries();
      if (!simFile.IsNull()) {
         fS = TFile::Open(simFile);
         if (fS && !fS->IsZombie()) {
            fT_S = (TTree *)fS->Get("cbmsim");
            fMcPts = new TClonesArray("AtMCPoint");
            fT_S->SetBranchAddress("AtTpcPoint", &fMcPts);
            fHasTruth = true;
         }
      }
      std::cout << "PRAExplorer: " << fNEvents << " events; truth=" << (fHasTruth ? "yes" : "no") << "\n";
      MakeGui();
      Goto(0);
   }

   void Next() { Goto(fCurrent + 1); }
   void Prev() { Goto(fCurrent - 1); }
   void GotoFromEntry() { if (fEvNum) Goto((int)fEvNum->GetNumberEntry()->GetIntNumber()); }
   void Goto(int e)
   {
      if (e < 0) e = 0;
      if (e >= fNEvents) e = fNEvents - 1;
      fCurrent = e;
      if (fEvNum) fEvNum->GetNumberEntry()->SetIntNumber(fCurrent);
      RunAndDraw();
   }
   void Run() { RunAndDraw(); } // re-run PRA on current event with current params

private:
   TGNumberEntry *addNum(TGCompositeFrame *p, const char *lab, double val, int digits = 4)
   {
      auto *row = new TGHorizontalFrame(p);
      auto *l = new TGLabel(row, lab);
      l->SetTextJustify(kTextLeft);
      row->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 1, 1));
      auto *ne = new TGNumberEntry(row, val, 7, -1, TGNumberFormat::kNESReal,
                                   TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
      row->AddFrame(ne, new TGLayoutHints(kLHintsRight, 2, 2, 1, 1));
      p->AddFrame(row, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      return ne;
   }
   TGCheckButton *addChk(TGCompositeFrame *p, const char *lab, bool on)
   {
      auto *cb = new TGCheckButton(p, lab);
      cb->SetState(on ? kButtonDown : kButtonUp);
      p->AddFrame(cb, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 1, 1));
      return cb;
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1500, 820);
      main->SetWindowName("PRA Explorer");
      auto *body = new TGHorizontalFrame(main);

      // ---- left: parameter panel ----
      auto *panel = new TGVerticalFrame(body, 240);
      panel->AddFrame(new TGLabel(panel, "--- TriplClust ---"), new TGLayoutHints(kLHintsTop, 2, 2, 3, 1));
      nT = addNum(panel, "t Tcluster", 8.0);
      nS = addNum(panel, "s Scluster", 0.3);
      nK = addNum(panel, "k Ktriplet", 19);
      nN = addNum(panel, "n Ntriplet", 2);
      nM = addNum(panel, "m Mcluster", 15);
      nR = addNum(panel, "r Rsmooth", 2.0);
      nA = addNum(panel, "a Atriplet", 0.03);
      panel->AddFrame(new TGLabel(panel, "--- clustering ---"), new TGLayoutHints(kLHintsTop, 2, 2, 4, 1));
      nCR = addNum(panel, "clusterRadius", 20);
      nCD = addNum(panel, "clusterDist", 15);
      nTC = addNum(panel, "arcTargetClus", 8);
      cArc = addChk(panel, "use ArcWalk", true);
      panel->AddFrame(new TGLabel(panel, "--- merge / charge ---"), new TGLayoutHints(kLHintsTop, 2, 2, 4, 1));
      cCircle = addChk(panel, "CircleMerge", true);
      nCRT = addNum(panel, " circ Rtol", 0.2);
      nCCT = addNum(panel, " circ Ctol mm", 25);
      cSelMerge = addChk(panel, "SelectAndMerge", false);
      cCharge = addChk(panel, "chargeFromCenter", true);
      nMH = addNum(panel, "minHits/track", 0);

      auto *bRun = new TGTextButton(panel, "  Run PRA  ");
      bRun->Connect("Clicked()", "PRAExplorer", this, "Run()");
      panel->AddFrame(bRun, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 8, 2));
      body->AddFrame(panel, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 2, 2, 2));

      // ---- right: nav bar + canvas ----
      auto *right = new TGVerticalFrame(body);
      auto *bar = new TGHorizontalFrame(right);
      auto *bPrev = new TGTextButton(bar, "<< Prev");
      bPrev->Connect("Clicked()", "PRAExplorer", this, "Prev()");
      bar->AddFrame(bPrev, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bNext = new TGTextButton(bar, "Next >>");
      bNext->Connect("Clicked()", "PRAExplorer", this, "Next()");
      bar->AddFrame(bNext, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      bar->AddFrame(new TGLabel(bar, "Event:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 4, 4));
      fEvNum = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0, 1e9);
      bar->AddFrame(fEvNum, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bGoto = new TGTextButton(bar, "Goto");
      bGoto->Connect("Clicked()", "PRAExplorer", this, "GotoFromEntry()");
      bar->AddFrame(bGoto, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      right->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *ec = new TRootEmbeddedCanvas("ec_pra", right, 1240, 760);
      right->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();

      body->AddFrame(right, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
      main->AddFrame(body, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   void RunAndDraw()
   {
      fT_D->GetEntry(fCurrent);
      if (fHasTruth) fT_S->GetEntry(fCurrent);
      fCanvas->cd();
      fCanvas->Clear();

      auto *ev = (AtEvent *)fEvArr->At(0);
      if (ev == nullptr) { fCanvas->Update(); return; }

      // configure the finder from the current GUI values
      AtPATTERN::AtTrackFinderTC tc;
      tc.SetTcluster(nT->GetNumber());
      tc.SetScluster(nS->GetNumber());
      tc.SetKtriplet((size_t)nK->GetNumber());
      tc.SetNtriplet((size_t)nN->GetNumber());
      tc.SetMcluster((size_t)nM->GetNumber());
      tc.SetRsmooth(nR->GetNumber());
      tc.SetAtriplet(nA->GetNumber());
      tc.SetClusterRadius(nCR->GetNumber());
      tc.SetClusterDistance(nCD->GetNumber());
      tc.SetUseArcWalk(cArc->IsOn());
      tc.SetTargetClusters((int)nTC->GetNumber());
      tc.SetUseSelectAndMerge(cSelMerge->IsOn());
      tc.SetUseCircleMerge(cCircle->IsOn());
      tc.SetCircleMergeTol(nCRT->GetNumber(), nCCT->GetNumber());
      tc.SetChargeFromCenter(cCharge->IsOn());
      tc.SetDiffusionParams(1e-5, 3e-4, 1.5, 0.08, 2.0); // PUMA-ish, so arc-walk cov is sane

      std::unique_ptr<AtPatternEvent> pe;
      try {
         pe = tc.FindTracks(*ev);
      } catch (...) {
         std::cerr << "PRA threw on event " << fCurrent << "\n";
      }

      // ---- draw: XY pad plane ----
      auto *fr = gPad->DrawFrame(-130, -130, 130, 130,
                                 Form("event %d  (edit params -> Run PRA);x [mm];y [mm]", fCurrent));
      (void)fr;
      auto ring = [](double r) { auto *e = new TEllipse(0, 0, r, r); e->SetFillStyle(0);
         e->SetLineColor(kGray + 1); e->SetLineStyle(2); e->Draw(); };
      ring(62.9); ring(121.1);

      // all PSA hits (grey)
      int nh = ev->GetNumHits();
      if (nh > 0) { auto *g = new TPolyMarker(nh);
         for (int i = 0; i < nh; ++i) { auto p = ev->GetHit(i).GetPosition(); g->SetPoint(i, p.X(), p.Y()); }
         g->SetMarkerColor(kGray + 2); g->SetMarkerStyle(1); g->Draw(); }

      // MC truth points (black open circles) if available
      if (fHasTruth && fMcPts->GetEntries() > 0) { int nm = fMcPts->GetEntries();
         auto *g = new TPolyMarker(nm);
         for (int k = 0; k < nm; ++k) { auto *mp = (AtMCPoint *)fMcPts->At(k); g->SetPoint(k, mp->GetX() * 10, mp->GetY() * 10); }
         g->SetMarkerColor(kBlack); g->SetMarkerStyle(24); g->SetMarkerSize(0.35); g->Draw(); }

      // tracks: hits coloured + fitted circle
      int cols[10] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                      kCyan + 2, kSpring - 6, kViolet + 1, kPink + 6, kAzure + 7};
      auto *info = new TPaveText(0.005, 0.005, 0.38, 0.30, "NDC");
      info->SetFillColor(0); info->SetTextAlign(12); info->SetTextSize(0.022);
      int nTrk = 0;
      if (pe) {
         auto &tc_arr = pe->GetTrackCand();
         nTrk = tc_arr.size();
         info->AddText(Form("tracks: %d", nTrk));
         int idx = 0;
         for (auto &tr : tc_arr) {
            int col = cols[idx % 10];
            auto &hits = tr.GetHitArray();
            if (!hits.empty()) { auto *g = new TPolyMarker(hits.size());
               for (size_t i = 0; i < hits.size(); ++i) { auto p = hits[i]->GetPosition(); g->SetPoint(i, p.X(), p.Y()); }
               g->SetMarkerColor(col); g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->Draw(); }
            double R = tr.GetGeoRadius();
            auto c = tr.GetGeoCenter();
            if (R > 0 && R < 1e5) { auto *e = new TEllipse(c.first, c.second, R, R);
               e->SetFillStyle(0); e->SetLineColor(col); e->SetLineWidth(2); e->Draw();
               double p = 0.299792458 * 4.0 * R;
               TString s = Form("#color[%d]{trk%d: R=%.0f  p=%.0f MeV", col, idx, R, p);
               if (fP0 > 0) s += Form("  (%+.0f%%)", 100 * (p - fP0) / fP0);
               s += Form("  %zu hits}", hits.size());
               info->AddText(s);
            } else {
               info->AddText(Form("#color[%d]{trk%d: R invalid (%zu hits)}", col, idx, hits.size()));
            }
            ++idx;
         }
      } else info->AddText("PRA: no result");
      info->Draw();
      std::cout << "event " << fCurrent << ": " << nh << " hits -> " << nTrk << " tracks\n";
      fCanvas->Update();
   }

   TFile *fD{nullptr}, *fS{nullptr};
   TTree *fT_D{nullptr}, *fT_S{nullptr};
   TClonesArray *fEvArr{nullptr}, *fMcPts{nullptr};
   int fCurrent{0}, fNEvents{0};
   bool fHasTruth{false};
   double fP0{0};
   TCanvas *fCanvas{nullptr};
   TGNumberEntry *fEvNum{nullptr};
   TGNumberEntry *nT, *nS, *nK, *nN, *nM, *nR, *nA, *nCR, *nCD, *nTC, *nCRT, *nCCT, *nMH;
   TGCheckButton *cCircle, *cSelMerge, *cCharge, *cArc;

   ClassDef(PRAExplorer, 0);
};

void pra_explorer(TString digiFile = "/mnt/f/puma_sweep/output_digi_pi200_primary.root", double p0 = 200.0,
                  TString simFile = "")
{
   gSystem->Load("libAtReconstruction.so");
   new PRAExplorer(digiFile, p0, simFile);
}
