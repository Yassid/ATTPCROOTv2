/// @file pid_compare_C15d.C
/// @brief Side-by-side PID plane, RAW vs GAIN MATCHED, with live binning and range controls.
///
///   root -l 'pid_compare_C15d.C()'
///   pid/open_pid_compare.sh          (detached, keeps the event loop alive)
///
/// WHY BOTH PANELS FROM ONE READ. The two planes must differ ONLY by the gain factor, so they are
/// filled from the same in-memory tracks in the same loop. Reading a "raw" plane from one file and
/// a "matched" plane from another invites them to differ by something else as well -- a different
/// run range, a different valid cut -- and then the comparison is not a comparison.
///
/// The gain is applied as sqrt(dE/dx) * sqrt(f), because the plane's x axis is sqrt(dE/dx) while
/// the factor f multiplies dE/dx itself. Getting that wrong shifts the matched band by sqrt(f)
/// instead of f and still looks plausible.
///
/// Controls: x and y bin counts, x and y maxima, log z, and a Brho slice projection that shows the
/// two band profiles overlaid -- which is where the matching is actually judged, since a 2D plane
/// at this statistics level hides a 10 % shift that a 1D profile makes obvious.

#include "gain_C15d.h"

class C15dPidCompare : public TObject {
public:
   C15dPidCompare(TString inDir, TString gainCsv, double xMax, double yMax, int nbx, int nby)
      : fInDir(inDir), fXmax(xMax), fYmax(yMax), fNbx(nbx), fNby(nby)
   {
      fGain = LoadGainTable_C15d(gainCsv.Data());
      if (fGain.empty())
         printf("\033[1;33mWARNING: no gain table at %s -- both panels will be identical\033[0m\n", gainCsv.Data());
      Load();
      MakeGui();
      Redraw();
   }

   /// read every <run>_pid.root once, keep the columns both panels need
   void Load()
   {
      TString ls = gSystem->GetFromPipe("ls -1 " + fInDir + "*_pid.root 2>/dev/null");
      TObjArray *files = ls.Tokenize("\n");
      long nAll = 0;
      for (int i = 0; i < files->GetEntries(); ++i) {
         TString fn = ((TObjString *)files->At(i))->GetString();
         TString base = gSystem->BaseName(fn.Data());
         // digits start at index 4 of "run_0013_pid.root". base(3,4) is "_001" and Atoi gives 0,
         // which silently maps every run onto run 0 -- that bug made the gain report claim the
         // whole set was a single run.
         const int run = TString(base(4, 4)).Atoi();
         TFile f(fn);
         auto *t = (TTree *)f.Get("pid");
         if (!t) continue;
         Double_t dedx = 0, brho = 0;
         Int_t valid = 1;
         t->SetBranchAddress("dEdx", &dedx);
         t->SetBranchAddress("brho", &brho);
         if (t->GetBranch("valid"))
            t->SetBranchAddress("valid", &valid);
         bool miss = false;
         const double f_ = GainFactor_C15d(fGain, run, miss);
         if (miss)
            ++fNmiss;
         for (Long64_t k = 0; k < t->GetEntries(); ++k) {
            t->GetEntry(k);
            ++nAll;
            if (!valid || dedx <= 0)
               continue;
            fSq.push_back(std::sqrt(dedx));
            fBr.push_back(brho);
            fF.push_back(f_);
            fRun.push_back(run);
         }
         ++fNruns;
      }
      printf("\033[1;32m  loaded %zu valid tracks of %ld, from %d runs (%d with no gain entry)\033[0m\n", fSq.size(),
             nAll, fNruns, fNmiss);
   }

   void Redraw()
   {
      const int nbx = (int)fEnbx->GetNumber(), nby = (int)fEnby->GetNumber();
      const double xhi = fExhi->GetNumber(), yhi = fEyhi->GetNumber();
      delete fHraw;
      delete fHcor;
      fHraw = new TH2F("hraw", "RAW  (no gain match);#sqrt{dE/dx};B#rho [T#upointm]", nbx, 0, xhi, nby, 0, yhi);
      fHcor = new TH2F("hcor", "GAIN MATCHED;#sqrt{dE/dx};B#rho [T#upointm]", nbx, 0, xhi, nby, 0, yhi);
      for (size_t i = 0; i < fSq.size(); ++i) {
         fHraw->Fill(fSq[i], fBr[i]);
         // the factor multiplies dE/dx, so on a sqrt(dE/dx) axis it enters as sqrt(f)
         fHcor->Fill(fSq[i] * std::sqrt(fF[i]), fBr[i]);
      }
      const bool logz = fLogZ->IsOn();
      fCanvas->Clear();
      fCanvas->Divide(2, 1);
      for (int p = 1; p <= 2; ++p) {
         fCanvas->cd(p);
         gPad->SetLogz(logz);
         gPad->SetRightMargin(0.13);
         (p == 1 ? fHraw : fHcor)->Draw("colz");
      }
      fCanvas->Update();
      fLabel->SetText(Form("  %zu tracks, %d runs   |   %d x %d bins, x < %.0f, y < %.2f", fSq.size(), fNruns, nbx,
                           nby, xhi, yhi));
   }

   /// the 1D test: one Brho slice, both profiles overlaid. A 2D plane hides a 10 % band shift.
   void Slice()
   {
      const double lo = fEblo->GetNumber(), hi = fEbhi->GetNumber();
      if (!(hi > lo)) {
         printf("Brho slice: hi must exceed lo\n");
         return;
      }
      const int nbx = (int)fEnbx->GetNumber();
      const double xhi = fExhi->GetNumber();
      auto *pr = new TH1F("prRaw", Form("B#rho %.2f-%.2f;#sqrt{dE/dx};tracks", lo, hi), nbx, 0, xhi);
      auto *pc = new TH1F("prCor", "", nbx, 0, xhi);
      for (size_t i = 0; i < fSq.size(); ++i) {
         if (fBr[i] < lo || fBr[i] >= hi)
            continue;
         pr->Fill(fSq[i]);
         pc->Fill(fSq[i] * std::sqrt(fF[i]));
      }
      auto *c = new TCanvas(Form("cslice%d", fSliceN++), "Brho slice", 900, 620);
      c->SetLogy(fLogZ->IsOn());
      pr->SetLineColor(kGray + 2);
      pr->SetLineWidth(2);
      pc->SetLineColor(kAzure + 2);
      pc->SetLineWidth(2);
      pr->Draw("hist");
      pc->Draw("hist same");
      auto *lg = new TLegend(0.58, 0.75, 0.88, 0.88);
      lg->AddEntry(pr, Form("raw  (RMS %.2f)", pr->GetRMS()), "l");
      lg->AddEntry(pc, Form("matched  (RMS %.2f)", pc->GetRMS()), "l");
      lg->Draw();
      c->Update();
      printf("  Brho %.2f-%.2f : raw RMS %.3f, matched RMS %.3f  -> %s\n", lo, hi, pr->GetRMS(), pc->GetRMS(),
             pc->GetRMS() < pr->GetRMS() ? "narrower" : "WIDER");
   }

   void SavePNG()
   {
      gSystem->mkdir("plots", kTRUE);
      fCanvas->SaveAs("plots/pid_compare_C15d.png");
      printf("  wrote plots/pid_compare_C15d.png\n");
   }
   void Quit() { gApplication->Terminate(0); }

private:
   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1500, 820);
      main->SetWindowName("C15d PID plane  --  raw vs gain matched (a2091 D2)");

      auto *bar = new TGHorizontalFrame(main);
      auto add = [&](const char *lab, TGNumberEntry *&e, double val, double mn, double mx) {
         bar->AddFrame(new TGLabel(bar, lab), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 2, 4, 4));
         e = new TGNumberEntry(bar, val, 6, -1, TGNumberFormat::kNESReal, TGNumberFormat::kNEANonNegative,
                               TGNumberFormat::kNELLimitMinMax, mn, mx);
         bar->AddFrame(e, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 4, 4));
      };
      add("x bins", fEnbx, fNbx, 20, 2000);
      add("y bins", fEnby, fNby, 20, 2000);
      add("x max", fExhi, fXmax, 5, 200);
      add("y max", fEyhi, fYmax, 0.2, 10);

      fLogZ = new TGCheckButton(bar, "log z");
      fLogZ->SetOn();
      bar->AddFrame(fLogZ, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 6, 4, 4));

      auto *bRe = new TGTextButton(bar, "Redraw");
      bRe->Connect("Clicked()", "C15dPidCompare", this, "Redraw()");
      bar->AddFrame(bRe, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));

      auto *bPng = new TGTextButton(bar, "Save PNG");
      bPng->Connect("Clicked()", "C15dPidCompare", this, "SavePNG()");
      bar->AddFrame(bPng, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));

      auto *bQ = new TGTextButton(bar, "Quit");
      bQ->Connect("Clicked()", "C15dPidCompare", this, "Quit()");
      bar->AddFrame(bQ, new TGLayoutHints(kLHintsLeft, 2, 6, 4, 4));
      main->AddFrame(bar, new TGLayoutHints(kLHintsExpandX));

      auto *bar2 = new TGHorizontalFrame(main);
      bar2->AddFrame(new TGLabel(bar2, "Brho slice"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 2, 4, 4));
      fEblo = new TGNumberEntry(bar2, 0.30, 6, -1, TGNumberFormat::kNESReal);
      bar2->AddFrame(fEblo, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 3, 4, 4));
      fEbhi = new TGNumberEntry(bar2, 0.40, 6, -1, TGNumberFormat::kNESReal);
      bar2->AddFrame(fEbhi, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 4, 4));
      auto *bSl = new TGTextButton(bar2, "Project both bands");
      bSl->Connect("Clicked()", "C15dPidCompare", this, "Slice()");
      bar2->AddFrame(bSl, new TGLayoutHints(kLHintsLeft, 4, 8, 4, 4));
      fLabel = new TGLabel(bar2, "");
      bar2->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
      main->AddFrame(bar2, new TGLayoutHints(kLHintsExpandX));

      auto *ec = new TRootEmbeddedCanvas("ec", main, 1480, 700);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();

      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   TString fInDir;
   double fXmax, fYmax;
   int fNbx, fNby, fNruns{0}, fNmiss{0}, fSliceN{0};
   std::map<int, double> fGain;
   std::vector<float> fSq, fBr, fF;
   std::vector<int> fRun;
   TH2F *fHraw{nullptr}, *fHcor{nullptr};
   TCanvas *fCanvas{nullptr};
   TGNumberEntry *fEnbx{nullptr}, *fEnby{nullptr}, *fExhi{nullptr}, *fEyhi{nullptr};
   TGNumberEntry *fEblo{nullptr}, *fEbhi{nullptr};
   TGCheckButton *fLogZ{nullptr};
   TGLabel *fLabel{nullptr};

   ClassDef(C15dPidCompare, 0)
};

void pid_compare_C15d(TString inDir = "/home/yassid/C15d_reco/", TString gainCsv = "gainmatch_C15d.csv",
                      double xMax = 85.0, double yMax = 2.5, int nbx = 340, int nby = 300)
{
   gStyle->SetOptStat(0);
   new C15dPidCompare(inDir, gainCsv, xMax, yMax, nbx, nby);
}
