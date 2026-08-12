/// @file gate_draw_pp.C
/// @brief Draw the a1975 (p,p) proton gate by hand, with the simulation telling you what it costs.
///
/// The point of this tool over a bare TCutG: when you close a polygon it immediately reports
///   - what fraction of TRUTH-MATCHED SIMULATED PROTONS the gate keeps  <- the efficiency
///   - how many DATA tracks it admits, against the current production gate
///   - how many of those also sit inside the deuteron gate               <- the contamination
/// Drawing on the data plane alone cannot show any of that, which is how the production gate ended
/// up keeping 61 % of protons, cutting 98 % of everything above 30 MeV, and carrying a notch near
/// sqrt(dE/dx) 8-11 that lands on the theta_cm 50-60 deg region.
///
/// The simulated proton band is drawn as CONTOURS over the data colz, so you can see where real
/// protons live while drawing on the distribution you actually have to cut.
///
///   root -l 'pid/gate_draw_pp.C'                         // defaults, writes pid/proton_hand.json
///   root -l 'pid/gate_draw_pp.C(106,140,"pid/my_gate.json")'
///
/// Buttons: [Draw new gate] -> click vertices, DOUBLE-CLICK to close -> [Evaluate] -> [Save JSON].
/// The sim points are cached to pid/plots/sim_proton_points.root on first run (~1 min); delete
/// that file to rebuild it after new simulation seeds.

// Explicit, rather than relying on the gSystem->Load in the constructor: cling parses the whole
// file before any of it runs, so without these the member accesses below are incomplete types.
#include "AtMCTrack.h"
#include "AtPIDEvent.h"
#include "AtSpyralPID.h"

#include <vector>

class PPGateDraw : public TObject {
public:
   PPGateDraw(int runLo, int runHi, TString outJson, TString dataDir, TString dataSuffix, TString simDir,
              TString simTags, TString curGate, TString deutGate, TString simCache, TString dataCache)
      : fOut(outJson), fSimCache(simCache), fDataCache(dataCache)
   {
      gSystem->Load("libAtReconstruction.so");
      gSystem->Load("libAtSimulationData.so");
      fH = new TH2F("hpp", "^{16}C(p,p) PID -- data colz, simulated protons in contour;"
                           "#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
      fS = new TH2F("hsim", "sim", 150, 0, 40, 125, 0, 1.2);

      // ---- data ----
      // A cache built by pid/make_data_points.C takes precedence. That is the only way to see the
      // plane at an fMinPoints other than 30: the AtPIDEvent persisted in the genfitter files was
      // written at the class default and cannot be re-cut, because the tracks it rejected have no
      // entry at all. Falls back to the persisted branch when no cache is given.
      long n = 0;
      if (fDataCache.Length() && !gSystem->AccessPathName(fDataCache)) {
         TFile *f = TFile::Open(fDataCache);
         auto *t = (TTree *)f->Get("pts");
         if (t) {
            float x, y; t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i); fDX.push_back(x); fDY.push_back(y); fH->Fill(x, y); ++n;
            }
            printf("  data: %ld PID points (from cache %s)\n", n, fDataCache.Data());
         }
         f->Close();
      }
      for (int r = runLo; n == 0 && r <= runHi; ++r) {
         TString fn = Form("%srun_%04d%s.root", dataDir.Data(), r, dataSuffix.Data());
         if (gSystem->AccessPathName(fn)) continue;
         TFile *f = TFile::Open(fn);
         auto *t = (TTree *)f->Get("cbmsim");
         if (!t) { f->Close(); continue; }
         t->SetBranchStatus("*", 0);
         t->SetBranchStatus("AtPIDEvent*", 1); // the tracking branch is large and unused here
         TClonesArray *pe = nullptr;
         t->SetBranchAddress("AtPIDEvent", &pe);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!pe || pe->GetEntries() == 0) continue;
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (!ev) continue;
            for (auto &s : ev->GetSpyral()) {
               if (!s.valid) continue;
               fH->Fill(s.sqrtdEdx, s.brho);
               fDX.push_back(s.sqrtdEdx); fDY.push_back(s.brho); ++n;
            }
         }
         f->Close();
      }
      if (!fDataCache.Length()) printf("  data: %ld PID points (runs %d-%d)\n", n, runLo, runHi);

      LoadSim(simDir, simTags);
      fCurN = Load(curGate, fCur, kRed + 1, 1);
      fDeuN = Load(deutGate, fDeu, kAzure + 2, 2);
      // how many data tracks the CURRENT production gate admits, for the comparison line
      fDataCur = 0;
      for (size_t i = 0; i < fDX.size(); ++i) if (In(fCurX, fCurY, fDX[i], fDY[i])) ++fDataCur;
      printf("  current gate admits %ld data tracks, keeps %.1f %% of simulated protons\n\n", fDataCur,
             fSX.empty() ? 0.0 : 100.0 * CountSim(fCurX, fCurY) / fSX.size());
      MakeGui();
      Redraw();
   }

   void DrawGate()
   {
      fCanvas->cd();
      printf("\n>>> click polygon vertices, DOUBLE-CLICK to close\n");
      TCutG *c = (TCutG *)fCanvas->WaitPrimitive("CUTG", "CutG");
      if (!c) { printf("no cut drawn\n"); return; }
      fNew = (TCutG *)c->Clone("newgate");
      fNew->SetLineColor(kGreen + 2);
      fNew->SetLineWidth(3);
      // The interactively drawn cut belongs to the pad. Deleting it without unlinking leaves a
      // dangling pointer in the primitive list, and the next Redraw walks that list -- which is
      // the "realloc(): invalid next size" heap corruption that lost a finished gate.
      if (fCanvas->GetListOfPrimitives()) fCanvas->GetListOfPrimitives()->Remove(c);
      delete c;
      // Save BEFORE anything else can crash. A polygon that exists only in memory is one segfault
      // away from being lost, and redrawing it by hand is not reproducible.
      Save();
      Evaluate();
      Redraw();
   }

   void Evaluate()
   {
      if (!fNew) { printf("draw a gate first\n"); return; }
      long sIn = 0;
      for (size_t i = 0; i < fSX.size(); ++i) if (fNew->IsInside(fSX[i], fSY[i])) ++sIn;
      long dIn = 0, dDeu = 0;
      for (size_t i = 0; i < fDX.size(); ++i)
         if (fNew->IsInside(fDX[i], fDY[i])) {
            ++dIn;
            if (fDeuN > 2 && In(fDeuX, fDeuY, fDX[i], fDY[i])) ++dDeu;
         }
      double eff = fSX.empty() ? 0 : 100.0 * sIn / fSX.size();
      double dRel = fDataCur ? 100.0 * (dIn - (double)fDataCur) / fDataCur : 0;
      TString msg = Form("  protons kept %.1f %% (%ld/%zu)   |   data %ld tracks (%+.1f %% vs current)"
                         "   |   deuteron overlap %.2f %%",
                         eff, sIn, fSX.size(), dIn, dRel, dIn ? 100.0 * dDeu / dIn : 0.0);
      fLabel->SetText(msg);
      printf("%s\n", msg.Data());
   }

   void Save()
   {
      if (!fNew) { printf("no new gate drawn -- nothing to save\n"); return; }
      FILE *f = fopen(fOut.Data(), "w");
      if (!f) { printf("cannot write %s\n", fOut.Data()); return; }
      fprintf(f, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n    \"vertices\": [\n",
              fName->GetText());
      int np = fNew->GetN();
      for (int i = 0; i < np; ++i) { double x, y; fNew->GetPoint(i, x, y);
         fprintf(f, "        [%.3f, %.3f]%s\n", x, y, (i == np - 1) ? "" : ","); }
      fprintf(f, "    ]\n}\n");
      fclose(f);
      printf("saved %d-vertex gate -> %s\n", np, fOut.Data());
      fLabel->SetText(Form("  saved %d vertices -> %s", np, fOut.Data()));
   }
   void SavePNG() { TString p = "pid/plots/gate_draw_pp.png"; fCanvas->SaveAs(p); printf("saved %s\n", p.Data()); }

   void Redraw()
   {
      fCanvas->cd(); fCanvas->SetLogz(); fCanvas->SetRightMargin(0.13);
      fH->Draw("colz");
      if (fS->GetEntries() > 0) {
         fS->SetContour(4);
         fS->SetLineColor(kGray + 2);
         fS->Draw("cont3 same");
      }
      if (fCur) fCur->Draw("L");
      if (fDeu) fDeu->Draw("L");
      if (fNew) fNew->Draw("L");
      auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.028);
      tx->SetTextColor(kRed + 1);    tx->DrawLatex(0.16, 0.88, "current production gate");
      tx->SetTextColor(kAzure + 2);  tx->DrawLatex(0.16, 0.85, "deuteron gate (avoid)");
      tx->SetTextColor(kGray + 2);   tx->DrawLatex(0.16, 0.82, "simulated protons (contour)");
      tx->SetTextColor(kGreen + 2);  tx->DrawLatex(0.16, 0.79, fNew ? "your gate" : "your gate: none yet");
      fCanvas->Modified(); fCanvas->Update(); gSystem->ProcessEvents();
   }

private:
   // truth-matched simulated protons, cached because the match needs the 225 MB sim files
   void LoadSim(TString simDir, TString simTags)
   {
      TString cache = fSimCache;
      if (!gSystem->AccessPathName(cache)) {
         TFile *f = TFile::Open(cache);
         auto *t = (TTree *)f->Get("pts");
         if (t) {
            float x, y; t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); fSX.push_back(x); fSY.push_back(y); fS->Fill(x, y); }
            printf("  sim: %zu truth-matched protons (from cache %s)\n", fSX.size(), cache.Data());
            f->Close();
            return;
         }
         f->Close();
      }
      printf("  building the simulated-proton cache (one minute, then it is reused)...\n");
      const double mp = 938.272;
      TObjArray *ta = simTags.Tokenize(",");
      for (int it = 0; it < ta->GetEntries(); ++it) {
         TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
         TString fs_ = simDir + "/" + tg + "_sim.root", ff_ = simDir + "/" + tg + "_genfitter_pp.root";
         if (gSystem->AccessPathName(fs_) || gSystem->AccessPathName(ff_)) continue;
         TFile *fs = TFile::Open(fs_), *ff = TFile::Open(ff_);
         TTree *ts = (TTree *)fs->Get("cbmsim"), *tf = (TTree *)ff->Get("cbmsim");
         if (!ts || !tf) { fs->Close(); ff->Close(); continue; }
         TClonesArray *mc = nullptr, *pe = nullptr;
         ts->SetBranchAddress("MCTrack", &mc);
         tf->SetBranchStatus("*", 0); tf->SetBranchStatus("AtPIDEvent*", 1);
         tf->SetBranchAddress("AtPIDEvent", &pe);
         for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
            ts->GetEntry(i); tf->GetEntry(i);
            double thT = -1;
            for (int k = 0; k < mc->GetEntriesFast(); ++k) {
               auto *p = (AtMCTrack *)mc->At(k);
               if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
               double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
               double pp = std::sqrt(px * px + py * py + pz * pz);
               if (pp <= 0) continue;
               thT = std::acos(pz / pp) * TMath::RadToDeg();
               break;
            }
            if (thT < 0 || !pe || pe->GetEntriesFast() == 0) continue;
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (!ev) continue;
            double bd = 1e9, bx = 0, by = 0; bool got = false;
            for (auto &sp : ev->GetSpyral()) {
               if (!sp.valid) continue;
               double d = std::fabs(180.0 - sp.polar * TMath::RadToDeg() - thT);
               if (d < bd) { bd = d; bx = sp.sqrtdEdx; by = sp.brho; got = true; }
            }
            if (got && bd < 10.0) { fSX.push_back(bx); fSY.push_back(by); fS->Fill(bx, by); }
         }
         fs->Close(); ff->Close();
      }
      gSystem->mkdir("pid/plots", kTRUE);
      TFile out(cache, "RECREATE");
      TTree t("pts", "truth-matched simulated protons");
      float x, y; t.Branch("x", &x); t.Branch("y", &y);
      for (size_t i = 0; i < fSX.size(); ++i) { x = fSX[i]; y = fSY[i]; t.Fill(); }
      t.Write(); out.Close();
      printf("  sim: %zu truth-matched protons (cached to %s)\n", fSX.size(), cache.Data());
   }

   int Load(TString fn, TPolyLine *&pl, int col, int sty)
   {
      std::vector<double> X, Y;
      if (!gSystem->AccessPathName(fn)) {
         std::ifstream in(fn.Data());
         std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
         size_t p = all.find("vertices");
         while (p != std::string::npos) {
            size_t a = all.find('[', p + 1); if (a == std::string::npos) break;
            size_t b = all.find(']', a);     if (b == std::string::npos) break;
            double x, y; char c;
            std::istringstream is(all.substr(a + 1, b - a - 1));
            if (is >> x >> c >> y) { X.push_back(x); Y.push_back(y); }
            p = b;
         }
      }
      if (X.size() > 2) {
         pl = new TPolyLine(X.size() + 1);
         for (size_t i = 0; i < X.size(); ++i) pl->SetPoint(i, X[i], Y[i]);
         pl->SetPoint(X.size(), X[0], Y[0]);
         pl->SetLineColor(col); pl->SetLineWidth(3); pl->SetLineStyle(sty);
      }
      if (fCurX.empty()) { fCurX = X; fCurY = Y; } else { fDeuX = X; fDeuY = Y; }
      return X.size();
   }
   static bool In(const std::vector<double> &X, const std::vector<double> &Y, double x, double y)
   {
      bool in = false; size_t n = X.size();
      if (n < 3) return false;
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((Y[i] > y) != (Y[j] > y)) && (x < (X[j] - X[i]) * (y - Y[i]) / (Y[j] - Y[i]) + X[i])) in = !in;
      return in;
   }
   long CountSim(const std::vector<double> &X, const std::vector<double> &Y)
   {
      long n = 0;
      for (size_t i = 0; i < fSX.size(); ++i) if (In(X, Y, fSX[i], fSY[i])) ++n;
      return n;
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1200, 950);
      main->SetWindowName("a1975 (p,p) proton gate");
      auto *bar = new TGHorizontalFrame(main);
      auto add = [&](const char *t, const char *slot) {
         auto *b = new TGTextButton(bar, t);
         b->Connect("Clicked()", "PPGateDraw", this, slot);
         bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      };
      add("  Draw new gate  ", "DrawGate()");
      add("  Evaluate  ", "Evaluate()");
      add("  Save JSON  ", "Save()");
      add("  Save PNG  ", "SavePNG()");
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fName = new TGTextEntry(bar, "proton_hand"); fName->Resize(150, 22);
      bar->AddFrame(fName, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  draw a polygon around the proton band; double-click to close, then Evaluate  ");
      main->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 6, 6, 2, 2));
      auto *ec = new TRootEmbeddedCanvas("ec_ppgate", main, 1180, 860);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      main->MapSubwindows(); main->Resize(main->GetDefaultSize()); main->MapWindow();
   }

   std::vector<float> fDX, fDY, fSX, fSY;
   std::vector<double> fCurX, fCurY, fDeuX, fDeuY;
   TH2F *fH{nullptr}, *fS{nullptr};
   TPolyLine *fCur{nullptr}, *fDeu{nullptr};
   TCutG *fNew{nullptr};
   TCanvas *fCanvas{nullptr};
   TGTextEntry *fName{nullptr};
   TGLabel *fLabel{nullptr};
   TString fOut;
   long fDataCur{0};
   int fCurN{0}, fDeuN{0};
   TString fSimCache, fDataCache;
   ClassDef(PPGateDraw, 0);
};

void gate_draw_pp(int runLo = 106, int runHi = 130, TString outJson = "pid/proton_hand.json",
                  TString dataDir = "/mnt/f/a1975/reco/", TString dataSuffix = "_genfitter_pp",
                  TString simDir = "/mnt/f/a1975_C16_pp_pid", TString simTags = "s2001,s2002,s2003,s2004,s2005,s2006",
                  TString curGate = "pid/proton_band.json", TString deutGate = "pid/deuteron_tight.json",
                  TString simCache = "pid/plots/sim_proton_points.root", TString dataCache = "")
{
   new PPGateDraw(runLo, runHi, outJson, dataDir, dataSuffix, simDir, simTags, curGate, deutGate, simCache,
                  dataCache);
}
