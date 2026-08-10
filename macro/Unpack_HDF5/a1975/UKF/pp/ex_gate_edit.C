/// @file ex_gate_edit.C
/// @brief Draw a cut by hand on the Ex vs theta_cm plane and save it, for selecting one state.
///
/// Same pattern as pp/pid_gate_edit.C -- TGMainFrame plus a plain 2D TCanvas, which works over
/// X11/WSLg; no Eve and no OpenGL, both of which crash here.
///
/// The plane is built with the correction tuned in the browser explorer, so what is drawn on is
/// what the analysis will use:
///
///     theta -> theta - (360/kcDenom)*(KE - kcPivot)   applied BEFORE Ex and theta_cm
///     Ebeam, kcDenom, kcPivot and the chi2 cut are all arguments, defaulting to the saved values
///
/// WHY A HAND-DRAWN POLYGON RATHER THAN A BAND. An automatic +-n sigma band around a fitted ridge
/// fails exactly where it matters: near theta_cm 50-65 the ground state is weak and the 1.766 MeV
/// 2+ is the tallest thing in the slice, so the fit walks onto the wrong peak. A polygon drawn by
/// eye follows the locus through those bins and stops where the locus genuinely stops.
///
/// The fitted ridge and its +-2.5 sigma envelope are drawn underneath as a guide, so the polygon
/// can follow them where they are trustworthy and depart from them where they are not.
///
///   root -l 'pp/ex_gate_edit.C("/path/pp_kin.root")'
///   root -l 'pp/ex_gate_edit.C("/path/pp_kin.root",188,4000,1.5,5,"gs")'
///
/// Workflow: [Draw cut] -> click vertices, DOUBLE-CLICK to close -> [Count inside] -> [Save].
/// Saves a JSON of the polygon AND a .root holding it as a TCutG, which the extraction step reads.

#include <vector>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <TROOT.h>
#include <TInterpreter.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TRootEmbeddedCanvas.h>
#include <iostream>

class ExGateEditor : public TObject {
public:
   ExGateEditor(TString cache, double Ebeam, double kcDenom, double kcPivot, double chi2Max, double icLo,
                double icHi, TString name, TString outStem)
      : fName0(name), fOut(outStem)
   {
      // Resolve the output directory HERE and keep it. gInterpreter->GetCurrentMacroName() is
      // only valid while a macro is executing; called later from a button handler it returns
      // null and constructing a TString from that segfaults -- which is exactly what happened
      // the first time Save() was pressed.
      {
         const char *mn = gInterpreter->GetCurrentMacroName();
         fDir = (mn && *mn) ? TString(gSystem->DirName(mn)) + "/" : TString("./");
      }
      const double slope = 360.0 / kcDenom;
      const double u = 931.49401;
      const double mb = 16.0147013 * u, mt = 1.00782503 * u, m3 = 1.00782503 * u, m4 = 16.0147013 * u;
      auto om = [](double x, double y, double z) {
         return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
      };

      fH = new TH2F("hex", TString::Format("^{16}C(p,p')  E_{x} vs #theta_{cm}   "
                                           "(E_{beam} = %.1f, #theta corr %.3f#circ/MeV about %.2f MeV);"
                                           "#theta_{cm} [deg];E_{x} [MeV]",
                                           Ebeam, slope, kcPivot),
                    180, 0, 180, 260, -3, 6);

      TFile *f = TFile::Open(cache);
      if (!f || f->IsZombie()) {
         std::cout << "\033[1;31mcannot open " << cache << "\033[0m\n";
         return;
      }
      TTree *t = (TTree *)f->Get("pk");
      float ke, th, vz, c2, ic;
      t->SetBranchAddress("ke", &ke);
      t->SetBranchAddress("theta", &th);
      t->SetBranchAddress("vz", &vz);
      t->SetBranchAddress("chi2ndf", &c2);
      t->SetBranchAddress("ic", &ic);
      long n = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!(c2 < chi2Max && ic > icLo && ic < icHi))
            continue;
         double thc = (th - slope * (ke - kcPivot)) * TMath::DegToRad();
         double E1 = Ebeam + mb, E3 = ke + m3, E4 = E1 + mt - E3;
         double s = mb * mb + mt * mt + 2 * mt * E1, uu = mt * mt + m3 * m3 - 2 * mt * E3;
         double a = (std::cos(thc) * om(s, mb * mb, mt * mt) * om(uu, mt * mt, m3 * m3) -
                     (s - mb * mb - mt * mt) * (mt * mt + m3 * m3 - uu)) /
                       (2 * mt * mt) +
                    s + uu - mt * mt;
         if (a < 0)
            continue;
         double m4x = std::sqrt(a), ex = m4x - m4;
         double tt = mt * mt + m4x * m4x - 2 * mt * E4;
         double arg = (s * s + s * (2 * tt - mb * mb - mt * mt - m3 * m3 - m4x * m4x) +
                       (mb * mb - mt * mt) * (m3 * m3 - m4x * m4x)) /
                      (om(s, mb * mb, mt * mt) * om(s, m3 * m3, m4x * m4x));
         if (arg < -1 || arg > 1)
            continue;
         double cm = (TMath::Pi() - std::acos(arg)) * TMath::RadToDeg();
         fH->Fill(cm, ex);
         fX.push_back(cm);
         fY.push_back(ex);
         ++n;
      }
      f->Close();
      std::cout << "ExGateEditor: " << n << " protons on the plane (Ebeam " << Ebeam << ", chi2 < " << chi2Max
                << ")\n";

      // the fitted ridge, as a guide only
      TString g = fDir + "plots/gs_cut_pp.root";
      if (!gSystem->AccessPathName(g)) {
         TFile *fg = TFile::Open(g);
         for (auto nm : {"gs_mu", "gs_lo", "gs_hi"}) {
            auto *gr = (TGraph *)fg->Get(nm);
            if (!gr)
               continue;
            auto *c = (TGraph *)gr->Clone();
            c->SetLineColor(kGray + 2);
            c->SetLineWidth(2);
            c->SetLineStyle(TString(nm) == "gs_mu" ? 1 : 2);
            fGuide.push_back(c);
         }
         fg->Close();
         std::cout << "  ridge guide loaded from plots/gs_cut_pp.root (grey)\n";
      }
      MakeGui();
      Redraw();
   }

   void Redraw()
   {
      fCanvas->cd();
      fCanvas->SetLogz();
      fCanvas->SetRightMargin(0.13);
      fH->Draw("colz");
      for (auto *g : fGuide)
         g->Draw("L same");
      if (fNew)
         fNew->Draw("L");
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.028);
      tx->SetTextColor(kGray + 2);
      tx->DrawLatex(0.16, 0.87, "fitted ridge #pm2.5#sigma (guide)");
      tx->SetTextColor(kGreen + 2);
      tx->DrawLatex(0.16, 0.83, fNew ? "drawn cut (will be saved)" : "drawn cut: none yet");
      fCanvas->Modified();
      fCanvas->Update();
      gSystem->ProcessEvents();
   }

   void DrawGate()
   {
      fCanvas->cd();
      std::cout << "\n>>> Click the polygon vertices around the state. DOUBLE-CLICK to close it.\n";
      TCutG *c = (TCutG *)fCanvas->WaitPrimitive("CUTG", "CutG");
      if (!c) {
         std::cout << "no cut drawn\n";
         return;
      }
      fNew = (TCutG *)c->Clone("excut");
      fNew->SetLineColor(kGreen + 2);
      fNew->SetLineWidth(3);
      delete c;
      CountInside();
      Redraw();
   }

   void CountInside()
   {
      if (!fNew) {
         std::cout << "draw a cut first\n";
         return;
      }
      long in = 0;
      for (size_t i = 0; i < fX.size(); ++i)
         if (fNew->IsInside(fX[i], fY[i]))
            ++in;
      fLabel->SetText(Form("  %ld / %zu protons inside (%.1f%%)  ", in, fX.size(),
                           100.0 * in / std::max((size_t)1, fX.size())));
      std::cout << "inside: " << in << " / " << fX.size() << "\n";
   }

   void Save()
   {
      if (!fNew) {
         std::cout << "no cut drawn -- nothing to save\n";
         return;
      }
      TString nm = fNameEntry->GetText();
      TString here = fDir;
      gSystem->mkdir(here + "plots", kTRUE);
      // JSON, so the cut is readable and diffable outside ROOT
      FILE *f = fopen(here + fOut + ".json", "w");
      if (!f)
         std::cout << "\033[1;31mcannot write " << here + fOut + ".json\033[0m\n";
      if (f) {
         fprintf(f, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"theta_cm\",\n    \"yaxis\": \"ex\",\n"
                    "    \"vertices\": [\n",
                 nm.Data());
         for (int i = 0; i < fNew->GetN(); ++i) {
            double x, y;
            fNew->GetPoint(i, x, y);
            fprintf(f, "        [%.3f, %.4f]%s\n", x, y, (i == fNew->GetN() - 1) ? "" : ",");
         }
         fprintf(f, "    ]\n}\n");
         fclose(f);
      }
      // and as a TCutG, which the extraction reads directly with IsInside()
      TFile fo(here + fOut + ".root", "RECREATE");
      if (fo.IsZombie()) {
         std::cout << "\033[1;31mcannot write " << here + fOut + ".root\033[0m\n";
         return;
      }
      fNew->SetName(nm);
      fNew->Write();
      fo.Close();
      std::cout << "saved " << fNew->GetN() << " vertices -> " << here + fOut << ".{json,root}\n";
      fLabel->SetText(Form("  saved %d vertices -> %s.{json,root}  ", fNew->GetN(), fOut.Data()));
   }

   void SavePNG()
   {
      TString p = fDir + "plots/ex_gate_edit.png";
      fCanvas->SaveAs(p);
      std::cout << "saved " << p << "\n";
   }

private:
   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1200, 900);
      main->SetWindowName("Ex vs theta_cm gate editor");
      auto *bar = new TGHorizontalFrame(main);
      auto *bD = new TGTextButton(bar, "  Draw cut  ");
      bD->Connect("Clicked()", "ExGateEditor", this, "DrawGate()");
      bar->AddFrame(bD, new TGLayoutHints(kLHintsLeft, 6, 4, 4, 4));
      auto *bC = new TGTextButton(bar, "  Count inside  ");
      bC->Connect("Clicked()", "ExGateEditor", this, "CountInside()");
      bar->AddFrame(bC, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bS = new TGTextButton(bar, "  Save  ");
      bS->Connect("Clicked()", "ExGateEditor", this, "Save()");
      bar->AddFrame(bS, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      auto *bP = new TGTextButton(bar, "  Save PNG  ");
      bP->Connect("Clicked()", "ExGateEditor", this, "SavePNG()");
      bar->AddFrame(bP, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fNameEntry = new TGTextEntry(bar, fName0.Data());
      fNameEntry->Resize(140, 22);
      bar->AddFrame(fNameEntry, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  Draw cut -> click vertices, double-click to close -> Count inside -> Save  ");
      main->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 6, 6, 2, 2));
      auto *ec = new TRootEmbeddedCanvas("ec_exedit", main, 1180, 820);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }
   std::vector<float> fX, fY;
   std::vector<TGraph *> fGuide;
   TH2F *fH{nullptr};
   TCutG *fNew{nullptr};
   TCanvas *fCanvas{nullptr};
   TGTextEntry *fNameEntry{nullptr};
   TGLabel *fLabel{nullptr};
   TString fName0, fOut, fDir;
   ClassDef(ExGateEditor, 0);
};

void ex_gate_edit(TString cache, Double_t Ebeam = 188.0, Double_t kcDenom = 4000.0, Double_t kcPivot = 1.5,
                  Double_t chi2Max = 5.0, TString name = "gs", TString outStem = "plots/ex_cut_gs",
                  Double_t icLo = 1000, Double_t icHi = 1350)
{
   new ExGateEditor(cache, Ebeam, kcDenom, kcPivot, chi2Max, icLo, icHi, name, outStem);
}
