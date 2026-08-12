/// @file gate_draw_sim.C
/// @brief Draw a proton gate on the SIMULATED PID plane. Simulation only -- no data anywhere.
///
/// Every point on the plot is a truth-matched simulated proton, so the only number that matters is
/// what fraction of them your polygon keeps. Evaluate reports exactly that.
///
///   root -l 'gate_draw_sim.C'                        // writes pid_gate_sim.json here
///   root -l 'gate_draw_sim.C("mygate.json")'
///
/// [Draw new gate] -> click vertices, DOUBLE-CLICK to close -> [Evaluate] -> [Save JSON].

#include "AtMCTrack.h"
#include "AtPIDEvent.h"
#include "AtSpyralPID.h"

#include <vector>

class SimGateDraw : public TObject {
public:
   SimGateDraw(TString outJson, TString simDir, TString tags, TString cacheFile, TString refGate,
               TString species)
      : fOut(outJson), fCache(cacheFile), fSpecies(species)
   {
      gSystem->Load("libAtReconstruction.so");
      gSystem->Load("libAtSimulationData.so");
      fH = new TH2F("hsim", "simulated " + fSpecies + ";#sqrt{dE/dx} [arb];B#rho [T#upointm]",
                    300, 0, 40, 250, 0, 2.4);
      Load(simDir, tags);
      LoadRef(refGate);
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
      long in = 0;
      for (size_t i = 0; i < fX.size(); ++i) if (fNew->IsInside(fX[i], fY[i])) ++in;
      TString m = Form("  keeps %ld / %zu simulated tracks  =  %.1f %%", in, fX.size(),
                       fX.empty() ? 0.0 : 100.0 * in / fX.size());
      fLabel->SetText(m);
      printf("%s\n", m.Data());
   }

   void Save()
   {
      if (!fNew) { printf("no gate drawn -- nothing to save\n"); return; }
      FILE *f = fopen(fOut.Data(), "w");
      if (!f) { printf("cannot write %s\n", fOut.Data()); return; }
      fprintf(f, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n    \"vertices\": [\n",
              fName->GetText());
      int np = fNew->GetN();
      for (int i = 0; i < np; ++i) { double x, y; fNew->GetPoint(i, x, y);
         fprintf(f, "        [%.3f, %.3f]%s\n", x, y, (i == np - 1) ? "" : ","); }
      fprintf(f, "    ]\n}\n");
      fclose(f);
      printf("saved %d vertices -> %s\n", np, fOut.Data());
      fLabel->SetText(Form("  saved %d vertices -> %s", np, fOut.Data()));
   }
   void SavePNG() { fCanvas->SaveAs("plots/gate_draw_sim.png"); printf("saved plots/gate_draw_sim.png\n"); }

   void Redraw()
   {
      fCanvas->cd(); fCanvas->SetLogz(); fCanvas->SetRightMargin(0.13);
      fH->Draw("colz");
      if (fRef) fRef->Draw("L");
      if (fNew) fNew->Draw("L");
      fCanvas->Modified(); fCanvas->Update(); gSystem->ProcessEvents();
   }

   /// Overlay a reference polygon (typically the hand-drawn DATA gate) so the sim gate can be
   /// traced against it. The two gates are drawn on different planes and are NOT required to be
   /// identical -- each has to select protons where its own plane puts them. How far apart they
   /// end up is the gain comparison, so the reference is shown, never enforced.
   void LoadRef(TString fn)
   {
      if (!fn.Length() || gSystem->AccessPathName(fn)) return;
      std::vector<double> X, Y;
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
      if (X.size() < 3) { printf("  reference gate %s: unusable (%zu vertices)\n", fn.Data(), X.size()); return; }
      fRef = new TPolyLine(X.size() + 1);
      for (size_t i = 0; i < X.size(); ++i) fRef->SetPoint(i, X[i], Y[i]);
      fRef->SetPoint(X.size(), X[0], Y[0]);
      fRef->SetLineColor(kRed + 1); fRef->SetLineWidth(3); fRef->SetLineStyle(2);
      // what the reference would keep HERE, i.e. if the data gate were used on the sim plane
      long in_ = 0;
      for (size_t i = 0; i < fX.size(); ++i) {
         bool ins = false;
         for (size_t k = 0, j = X.size() - 1; k < X.size(); j = k++)
            if (((Y[k] > fY[i]) != (Y[j] > fY[i])) &&
                (fX[i] < (X[j] - X[k]) * (fY[i] - Y[k]) / (Y[j] - Y[k]) + X[k]))
               ins = !ins;
         if (ins) ++in_;
      }
      printf("  reference gate %s (%zu vertices) keeps %.1f %% of simulated protons\n", fn.Data(), X.size(),
             fX.empty() ? 0.0 : 100.0 * in_ / fX.size());
   }

private:
   void Load(TString simDir, TString tags)
   {
      TString cache = fCache;
      if (!gSystem->AccessPathName(cache)) {
         TFile *f = TFile::Open(cache);
         auto *t = (TTree *)f->Get("pts");
         if (t) {
            float x, y; t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); fX.push_back(x); fY.push_back(y); fH->Fill(x, y); }
            printf("  %zu simulated protons (cache)\n", fX.size());
            f->Close();
            return;
         }
         f->Close();
      }
      const double mp = 938.272;
      TObjArray *ta = tags.Tokenize(",");
      for (int it = 0; it < ta->GetEntries(); ++it) {
         TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
         TString a = simDir + "/" + tg + "_sim.root", b = simDir + "/" + tg + "_genfitter_pp.root";
         if (gSystem->AccessPathName(a) || gSystem->AccessPathName(b)) continue;
         TFile *fs = TFile::Open(a), *ff = TFile::Open(b);
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
            if (got && bd < 10.0) { fX.push_back(bx); fY.push_back(by); fH->Fill(bx, by); }
         }
         fs->Close(); ff->Close();
      }
      gSystem->mkdir("plots", kTRUE);
      TFile out(cache, "RECREATE");
      TTree t("pts", "truth-matched simulated protons");
      float x, y; t.Branch("x", &x); t.Branch("y", &y);
      for (size_t i = 0; i < fX.size(); ++i) { x = fX[i]; y = fY[i]; t.Fill(); }
      t.Write(); out.Close();
      printf("  %zu simulated protons (cached)\n", fX.size());
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1150, 900);
      main->SetWindowName("simulated proton gate");
      auto *bar = new TGHorizontalFrame(main);
      auto add = [&](const char *t, const char *slot) {
         auto *b = new TGTextButton(bar, t);
         b->Connect("Clicked()", "SimGateDraw", this, slot);
         bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      };
      add("  Draw new gate  ", "DrawGate()");
      add("  Evaluate  ", "Evaluate()");
      add("  Save JSON  ", "Save()");
      add("  Save PNG  ", "SavePNG()");
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fName = new TGTextEntry(bar, "proton_sim_hand"); fName->Resize(150, 22);
      bar->AddFrame(fName, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  draw a polygon around the band; double-click to close, then Evaluate  ");
      main->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 6, 6, 2, 2));
      auto *ec = new TRootEmbeddedCanvas("ec_simgate", main, 1130, 820);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      main->MapSubwindows(); main->Resize(main->GetDefaultSize()); main->MapWindow();
   }

   std::vector<float> fX, fY;
   TH2F *fH{nullptr};
   TCutG *fNew{nullptr};
   TPolyLine *fRef{nullptr};
   TString fCache, fSpecies;
   TCanvas *fCanvas{nullptr};
   TGTextEntry *fName{nullptr};
   TGLabel *fLabel{nullptr};
   TString fOut;
   ClassDef(SimGateDraw, 0);
};

/// @param cacheFile  truth-matched proton points. Build one at a chosen fMinPoints with
///                   make_sim_points.C; the default was written at the class default of 30.
/// @param refGate    optional polygon drawn on top, e.g. the hand-drawn DATA gate, so the two can
///                   be kept close. Shown for reference only -- the two planes are not required to
///                   coincide, and how far apart they sit is the gain comparison.
void gate_draw_sim(TString outJson = "pid_gate_sim.json", TString simDir = "/mnt/f/a1975_C16_pp_pid",
                   TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                   TString cacheFile = "plots/sim_proton_points.root", TString refGate = "",
                   TString species = "protons  ^{16}C(p,p')")
{
   new SimGateDraw(outJson, simDir, tags, cacheFile, refGate, species);
}
