/// @file kin_gui_Be10.C
/// @brief Interactive browser for the 10Be(t,p)12Be proton kinematics: KE_p vs theta_lab, with
///        the analytic level curves over it and a live readout of what any point on the plot means.
///
///   root -l 'kin_gui_Be10.C("/mnt/f/Be10_tp")'
///
/// Launch it DETACHED, or ROOT reads EOF on stdin and exits the moment it is backgrounded:
///   setsid bash -c 'tail -f /dev/null | root -l "kin_gui_Be10.C(\"/mnt/f/Be10_tp\")"' &
///
/// WHAT IT SHOWS. One big canvas instead of a sixth of a page, and the histogram is the point:
///   - the generated truth, the accepted truth, or the reconstructed sample, per configuration
///   - any subset of the four levels
///   - the exact two-body curves on top, so the difference between a level and the band it
///     produces is visible rather than inferred
///   - iso-theta_cm lines, because theta_cm is the variable the physics is quoted in and it runs
///     BACKWARDS across this plot (small theta_cm is large theta_lab)
///
/// THE READOUT IS THE USEFUL PART. Moving the cursor inverts the two-body kinematics at that
/// (theta_lab, KE) and reports the excitation energy and theta_cm it corresponds to. That turns
/// the plot from a picture into something you can interrogate: how far in KE is the 0+_2 from the
/// 2+ at this angle, how much Ex does 100 keV of proton energy buy here, where does a given level
/// stop being separable.
///
/// Buttons open the E_x projection and the median-KE profile for whatever is currently selected.

#include "RQ_OBJECT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGButton.h"
#include "TGComboBox.h"
#include "TGFrame.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGStatusBar.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRootEmbeddedCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLatex.h"

#include <cmath>
#include <cstdio>
#include <vector>

// ---- 10Be(t,p)12Be kinematics -------------------------------------------------------------
static const double GU = 931.49401;
static const double GM1 = 10.0135341 * GU, GM2 = 3.0160493 * GU, GM3 = 1.007825 * GU, GM4 = 12.0269221 * GU;
static const double GEB = 112.20;
static const int GNL = 4;
static const double GLEX[GNL] = {0.0, 2.109, 2.251, 2.715};
static const char *GLVN[GNL] = {"gs", "ex2109", "ex2251", "ex2715"};
static const char *GLJP[GNL] = {"0^{+} g.s.", "2^{+} 2.109", "0^{+}_{2} 2.251", "1^{-} 2.715"};
static const int GLCOL[GNL] = {kBlue, kGreen + 2, kMagenta, kOrange + 7};

static double g_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
/// (theta_lab, KE) -> {Ex, theta_cm in the DWBA convention}
static void g_inv(double th, double Ke, double &Ex, double &cm)
{
   Ex = NAN;
   cm = NAN;
   double Et1 = GEB + GM1, Et3 = Ke + GM3;
   double s = GM1 * GM1 + GM2 * GM2 + 2 * GM2 * Et1;
   double uu = GM2 * GM2 + GM3 * GM3 - 2 * GM2 * Et3;
   double a = (std::cos(th) * g_om2(s, GM1 * GM1, GM2 * GM2) * g_om2(uu, GM2 * GM2, GM3 * GM3) -
               (s - GM1 * GM1 - GM2 * GM2) * (GM2 * GM2 + GM3 * GM3 - uu)) / (2 * GM2 * GM2) + s + uu - GM2 * GM2;
   if (a <= 0) return;
   double m4x = std::sqrt(a);
   Ex = m4x - GM4;
   double Et4 = Et1 + GM2 - Et3;
   double t = GM2 * GM2 + m4x * m4x - 2 * GM2 * Et4;
   double c = (s * s + s * (2 * t - GM1 * GM1 - GM2 * GM2 - GM3 * GM3 - m4x * m4x) +
               (GM1 * GM1 - GM2 * GM2) * (GM3 * GM3 - m4x * m4x)) /
              (g_om2(s, GM1 * GM1, GM2 * GM2) * g_om2(s, GM3 * GM3, m4x * m4x));
   if (c >= -1 && c <= 1) cm = (TMath::Pi() - std::acos(c)) * TMath::RadToDeg();
}
/// theta_cm (DWBA) -> (theta_lab, KE) for residual mass m4
static bool g_fwd(double m4, double thcmA, double &thlab, double &Ke)
{
   const double thcm = TMath::Pi() - thcmA;
   double E1 = GEB + GM1, s = GM1 * GM1 + GM2 * GM2 + 2 * GM2 * E1, rs = std::sqrt(s);
   if (rs < GM3 + m4) return false;
   double pcm = g_om2(s, GM3 * GM3, m4 * m4) / (2 * rs), Ecm3 = std::sqrt(pcm * pcm + GM3 * GM3);
   double plab = std::sqrt(E1 * E1 - GM1 * GM1), beta = plab / (E1 + GM2), gam = 1 / std::sqrt(1 - beta * beta);
   Ke = gam * (Ecm3 + beta * pcm * std::cos(thcm)) - GM3;
   thlab = std::atan2(pcm * std::sin(thcm), gam * (pcm * std::cos(thcm) + beta * Ecm3));
   return true;
}

class KinGui {
   RQ_OBJECT("KinGui")
private:
   TGMainFrame *fMain{nullptr};
   TRootEmbeddedCanvas *fEc{nullptr};
   TGStatusBar *fBar{nullptr};
   TGComboBox *fCfg{nullptr}, *fSample{nullptr};
   TGCheckButton *fLev[GNL]{}, *fCurves{nullptr}, *fCmLines{nullptr}, *fLogZ{nullptr};
   TGNumberEntry *fCmLo{nullptr}, *fCmHi{nullptr};
   TString fRoot, fOut;
   std::vector<TString> fCfgs;
   TH2D *fH{nullptr};
   std::vector<TObject *> fOverlay;

public:
   KinGui(const TGWindow *p, TString root, TString outDir);
   virtual ~KinGui();
   void Redraw();
   void OnMove(Int_t ev, Int_t px, Int_t py, TObject *sel);
   void ProjectEx();
   void ProfileKE();
   void SavePng();
   void Fill(TH2D *h, TH1D *hex);
};

KinGui::KinGui(const TGWindow *p, TString root, TString outDir) : fRoot(root), fOut(outDir)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   fMain = new TGMainFrame(p, 1500, 900);
   fMain->SetWindowName("10Be(t,p)12Be  --  KE_p vs theta_lab");
   auto *hf = new TGHorizontalFrame(fMain, 1500, 860);

   // ---- controls ----
   auto *ctl = new TGVerticalFrame(hf, 250, 860);
   ctl->AddFrame(new TGLabel(ctl, "configuration"), new TGLayoutHints(kLHintsLeft, 4, 4, 8, 2));
   fCfg = new TGComboBox(ctl, 100);
   for (const char *c : {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"}) {
      // only offer a configuration whose exres files actually exist, so a half-finished campaign
      // cannot be browsed as if it were complete
      TString t = gSystem->GetFromPipe(
         TString::Format("ls %s/%s/exres_gs_s*_%s.root 2>/dev/null | head -1", fRoot.Data(), c, c));
      if (t.Strip(TString::kBoth).IsNull()) continue;
      fCfgs.push_back(c);
      fCfg->AddEntry(c, (int)fCfgs.size() - 1);
   }
   if (fCfgs.empty()) { fCfgs.push_back("b285_attpc"); fCfg->AddEntry("b285_attpc", 0); }
   fCfg->Select(0);
   fCfg->Resize(230, 22);
   ctl->AddFrame(fCfg, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 6));

   ctl->AddFrame(new TGLabel(ctl, "sample"), new TGLayoutHints(kLHintsLeft, 4, 4, 4, 2));
   fSample = new TGComboBox(ctl, 101);
   fSample->AddEntry("generated truth", 0);
   fSample->AddEntry("accepted truth", 1);
   fSample->AddEntry("reconstructed", 2);
   fSample->Select(0);
   fSample->Resize(230, 22);
   ctl->AddFrame(fSample, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 8));

   ctl->AddFrame(new TGLabel(ctl, "levels"), new TGLayoutHints(kLHintsLeft, 4, 4, 4, 2));
   const char *lbl[GNL] = {"0+  g.s.", "2+  2.109", "0+_2 2.251", "1-  2.715"};
   for (int l = 0; l < GNL; ++l) {
      fLev[l] = new TGCheckButton(ctl, lbl[l], 200 + l);
      fLev[l]->SetOn(kTRUE);
      ctl->AddFrame(fLev[l], new TGLayoutHints(kLHintsLeft, 12, 4, 1, 1));
   }

   ctl->AddFrame(new TGLabel(ctl, "theta_cm window [deg]"), new TGLayoutHints(kLHintsLeft, 4, 4, 10, 2));
   auto *cmf = new TGHorizontalFrame(ctl);
   fCmLo = new TGNumberEntry(cmf, 0, 5, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative);
   fCmHi = new TGNumberEntry(cmf, 180, 5, -1, TGNumberFormat::kNESRealOne, TGNumberFormat::kNEANonNegative);
   cmf->AddFrame(fCmLo, new TGLayoutHints(kLHintsLeft, 8, 4, 0, 0));
   cmf->AddFrame(fCmHi, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 0));
   ctl->AddFrame(cmf, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 6));

   fCurves = new TGCheckButton(ctl, "analytic level curves", 300);
   fCurves->SetOn(kTRUE);
   ctl->AddFrame(fCurves, new TGLayoutHints(kLHintsLeft, 12, 4, 6, 1));
   fCmLines = new TGCheckButton(ctl, "iso-#theta_{cm} markers", 301);
   fCmLines->SetOn(kTRUE);
   ctl->AddFrame(fCmLines, new TGLayoutHints(kLHintsLeft, 12, 4, 1, 1));
   fLogZ = new TGCheckButton(ctl, "log z", 302);
   fLogZ->SetOn(kTRUE);
   ctl->AddFrame(fLogZ, new TGLayoutHints(kLHintsLeft, 12, 4, 1, 8));

   auto *bRe = new TGTextButton(ctl, "&Redraw", 400);
   ctl->AddFrame(bRe, new TGLayoutHints(kLHintsExpandX, 8, 8, 4, 3));
   auto *bEx = new TGTextButton(ctl, "E_x projection", 401);
   ctl->AddFrame(bEx, new TGLayoutHints(kLHintsExpandX, 8, 8, 3, 3));
   auto *bPr = new TGTextButton(ctl, "median KE profile", 402);
   ctl->AddFrame(bPr, new TGLayoutHints(kLHintsExpandX, 8, 8, 3, 3));
   auto *bSv = new TGTextButton(ctl, "save PNG", 403);
   ctl->AddFrame(bSv, new TGLayoutHints(kLHintsExpandX, 8, 8, 3, 10));

   auto *help = new TGLabel(ctl, "move the cursor over the plot:\nit inverts the kinematics there\nand reports E_x and theta_cm");
   ctl->AddFrame(help, new TGLayoutHints(kLHintsLeft, 6, 4, 4, 4));

   hf->AddFrame(ctl, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 2, 2, 2));

   fEc = new TRootEmbeddedCanvas("ec", hf, 1230, 850);
   hf->AddFrame(fEc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
   fMain->AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   int parts[3] = {40, 30, 30};
   fBar = new TGStatusBar(fMain, 1500, 22);
   fBar->SetParts(parts, 3);
   fMain->AddFrame(fBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));

   bRe->Connect("Clicked()", "KinGui", this, "Redraw()");
   bEx->Connect("Clicked()", "KinGui", this, "ProjectEx()");
   bPr->Connect("Clicked()", "KinGui", this, "ProfileKE()");
   bSv->Connect("Clicked()", "KinGui", this, "SavePng()");
   fCfg->Connect("Selected(Int_t)", "KinGui", this, "Redraw()");
   fSample->Connect("Selected(Int_t)", "KinGui", this, "Redraw()");
   for (int l = 0; l < GNL; ++l) fLev[l]->Connect("Clicked()", "KinGui", this, "Redraw()");
   fCurves->Connect("Clicked()", "KinGui", this, "Redraw()");
   fCmLines->Connect("Clicked()", "KinGui", this, "Redraw()");
   fLogZ->Connect("Clicked()", "KinGui", this, "Redraw()");
   fEc->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "KinGui", this,
                            "OnMove(Int_t,Int_t,Int_t,TObject*)");

   fMain->MapSubwindows();
   fMain->Resize(fMain->GetDefaultSize());
   fMain->MapWindow();
   Redraw();
}

KinGui::~KinGui()
{
   fMain->Cleanup();
}

/// Fill the 2D (and, if asked, the E_x spectrum) for one selection. A FREE function, not a member,
/// so the batch renderer below draws from exactly the same code the GUI does -- a figure that came
/// from a second copy of this loop would be a figure of something slightly different.
///   samp 0 = generated truth (from the cache), 1 = accepted truth, 2 = reconstructed
static void kin_fill(const TString &root, const TString &outDir, const TString &cfg, int samp, const bool *levOn,
                     double cl, double ch, TH2D *h, TH1D *hex)
{
   for (int l = 0; l < GNL; ++l) {
      if (!levOn[l]) continue;
      if (samp == 0) {
         TFile *f = TFile::Open(outDir + "/kin_truth_Be10.root");
         if (!f || f->IsZombie()) continue;
         TTree *t = (TTree *)f->Get(Form("truth_%s", GLVN[l]));
         if (!t) { f->Close(); continue; }
         double th, ke;
         t->SetBranchAddress("th", &th);
         t->SetBranchAddress("ke", &ke);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double ex, cm;
            g_inv(th * TMath::DegToRad(), ke, ex, cm);
            if (!(cm >= cl && cm < ch)) continue;
            if (h) h->Fill(th, ke);
            if (hex) hex->Fill(ex);
         }
         f->Close();
      } else {
         TString fn = gSystem->GetFromPipe(TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1",
                                                           root.Data(), cfg.Data(), GLVN[l], cfg.Data()));
         fn = fn.Strip(TString::kBoth);
         if (fn.IsNull()) continue;
         TFile *f = TFile::Open(fn);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (!t) { if (f) f->Close(); continue; }
         double thT, keT, thR, keR, cmT, exR;
         t->SetBranchAddress("thTrue", &thT);
         t->SetBranchAddress("keTrue", &keT);
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keReco", &keR);
         t->SetBranchAddress("cmTrue", &cmT);
         t->SetBranchAddress("exReco", &exR);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!(cmT >= cl && cmT < ch)) continue;
            if (samp == 1) { if (h) h->Fill(thT, keT); if (hex) hex->Fill(GLEX[l]); }
            else           { if (h) h->Fill(thR, keR); if (hex) hex->Fill(exR); }
         }
         f->Close();
      }
   }
}

void KinGui::Fill(TH2D *h, TH1D *hex)
{
   bool on[GNL];
   for (int l = 0; l < GNL; ++l) on[l] = fLev[l]->IsOn();
   kin_fill(fRoot, fOut, fCfgs[fCfg->GetSelected()], fSample->GetSelected(), on, fCmLo->GetNumber(),
            fCmHi->GetNumber(), h, hex);
}

void KinGui::Redraw()
{
   TCanvas *c = fEc->GetCanvas();
   c->Clear();
   c->SetRightMargin(0.13);
   c->SetLeftMargin(0.10);
   c->SetTopMargin(0.08);
   if (fH) { delete fH; fH = nullptr; }
   for (auto *o : fOverlay) delete o;
   fOverlay.clear();

   TString cfg = fCfgs[fCfg->GetSelected()];
   const char *sn[3] = {"generated truth", "accepted truth", "reconstructed"};
   int samp = fSample->GetSelected();

   fH = new TH2D("hKin", Form("^{10}Be(t,p)^{12}Be   %s   %s;#theta_{lab} [deg];KE_{p} [MeV]",
                              samp == 0 ? "(field-independent)" : cfg.Data(), sn[samp]),
                 360, 0, 180, 275, 0, 55);
   fH->SetDirectory(nullptr);
   Fill(fH, nullptr);
   c->SetLogz(fLogZ->IsOn() ? 1 : 0);
   fH->Draw("colz");

   if (fCurves->IsOn()) {
      auto *leg = new TLegend(0.55, 0.66, 0.86, 0.90);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.028);
      for (int l = 0; l < GNL; ++l) {
         if (!fLev[l]->IsOn()) continue;
         auto *g = new TGraph();
         int n = 0;
         for (double a = 0.5; a <= 179.5; a += 0.25) {
            double th, ke;
            if (!g_fwd(GM4 + GLEX[l], a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
            g->SetPoint(n++, th * TMath::RadToDeg(), ke);
         }
         g->SetLineColor(GLCOL[l]);
         g->SetLineWidth(2);
         g->Draw("L SAME");
         leg->AddEntry(g, GLJP[l], "l");
         fOverlay.push_back(g);
      }
      leg->Draw();
      fOverlay.push_back(leg);
   }
   if (fCmLines->IsOn()) {
      // theta_cm runs BACKWARDS across this plot, which is the single most confusing thing about
      // inverse kinematics, so mark it explicitly rather than leaving it to be remembered
      auto *g = new TGraph();
      auto *tx = new TLatex();
      tx->SetTextSize(0.022);
      tx->SetTextColor(kGray + 3);
      int n = 0;
      for (double a : {5., 10., 20., 30., 45., 60., 90., 120., 150.}) {
         double th, ke;
         if (!g_fwd(GM4, a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
         g->SetPoint(n++, th * TMath::RadToDeg(), ke);
         tx->DrawLatex(th * TMath::RadToDeg() + 1.5, ke + 1.0, Form("%.0f#circ", a));
      }
      g->SetMarkerStyle(29);
      g->SetMarkerSize(1.4);
      g->SetMarkerColor(kGray + 3);
      g->Draw("P SAME");
      fOverlay.push_back(g);
      fOverlay.push_back(tx);
      auto *lb = new TLatex();
      lb->SetNDC();
      lb->SetTextSize(0.024);
      lb->SetTextColor(kGray + 3);
      // INSIDE the frame, lower left: that corner is always empty (a forward proton cannot be
      // slow) and NDC 0.93 sat on top of the pad title.
      lb->DrawLatex(0.14, 0.16, "stars: #theta_{cm} (DWBA) -- runs BACKWARDS across this plot");
      fOverlay.push_back(lb);
   }
   c->Modified();
   c->Update();
   fBar->SetText(Form("%s | %s | %.0f entries", cfg.Data(), sn[samp], fH->GetEntries()), 0);
}

void KinGui::OnMove(Int_t, Int_t px, Int_t py, TObject *)
{
   TCanvas *c = fEc->GetCanvas();
   double x = c->AbsPixeltoX(px), y = c->AbsPixeltoY(py);
   double th = c->PadtoX(x), ke = c->PadtoY(y);
   if (th < 0 || th > 180 || ke < 0 || ke > 55) return;
   double ex, cm;
   g_inv(th * TMath::DegToRad(), ke, ex, cm);
   fBar->SetText(Form("theta_lab %.2f deg   KE_p %.3f MeV", th, ke), 1);
   if (std::isnan(ex)) fBar->SetText("outside the kinematic locus", 2);
   else fBar->SetText(Form("=>  E_x %+.3f MeV   theta_cm %.1f deg", ex, std::isnan(cm) ? 0.0 : cm), 2);
}

void KinGui::ProjectEx()
{
   auto *c = new TCanvas(Form("cEx%d", (int)gRandom->Integer(100000)), "Ex projection", 1000, 650);
   auto *h = new TH1D("hExProj", ";E_{x}(^{12}Be) [MeV];counts", 300, -1.5, 4.5);
   Fill(nullptr, h);
   h->SetLineColor(kBlack);
   h->Draw();
   for (int l = 0; l < GNL; ++l) {
      if (!fLev[l]->IsOn()) continue;
      auto *ln = new TLine(GLEX[l], 0, GLEX[l], h->GetMaximum() * 1.05);
      ln->SetLineColor(GLCOL[l]);
      ln->SetLineStyle(2);
      ln->Draw();
   }
   c->Modified();
   c->Update();
}

void KinGui::ProfileKE()
{
   auto *c = new TCanvas(Form("cPr%d", (int)gRandom->Integer(100000)), "KE profile", 1000, 650);
   auto *h = new TH2D("hPr", ";#theta_{lab} [deg];KE_{p} [MeV]", 180, 0, 180, 275, 0, 55);
   h->SetDirectory(nullptr);
   Fill(h, nullptr);
   auto *g = new TGraph();
   int n = 0;
   for (int b = 1; b <= h->GetNbinsX(); ++b) {
      TH1D *p = h->ProjectionY("_p", b, b);
      if (p->GetEntries() < 10) { delete p; continue; }
      double q = 0.5, med;
      p->GetQuantiles(1, &med, &q);
      g->SetPoint(n++, h->GetXaxis()->GetBinCenter(b), med);
      delete p;
   }
   h->Draw("colz");
   g->SetLineColor(kRed);
   g->SetLineWidth(3);
   g->Draw("L SAME");
   c->Modified();
   c->Update();
}

void KinGui::SavePng()
{
   TString f = TString::Format("%s/kin_gui_%s_s%d.png", fOut.Data(), fCfgs[fCfg->GetSelected()].Data(),
                               fSample->GetSelected());
   fEc->GetCanvas()->SaveAs(f);
   fBar->SetText(Form("wrote %s", f.Data()), 0);
}

void kin_gui_Be10(TString root = "/mnt/f/Be10_tp", TString outDir = "plots")
{
   if (gSystem->AccessPathName(outDir + "/kin_truth_Be10.root"))
      printf("\033[1;31mno %s/kin_truth_Be10.root -- run make_kin_points_Be10.C first; the\n"
             "'generated truth' option will be empty until you do.\033[0m\n", outDir.Data());
   new KinGui(gClient->GetRoot(), root, outDir);
}

/// A single big publication-quality panel of the same plot, without the GUI -- for a figure, and
/// for checking what the GUI is showing without a display.
///   root -b -q 'kin_gui_Be10.C+' ... or simply:
///   root -b -q 'kin_png_Be10.C("/mnt/f/Be10_tp","b285_attpc",0)'
void kin_png_Be10(TString root = "/mnt/f/Be10_tp", TString cfg = "b285_attpc", int samp = 0,
                  TString outDir = "plots", double cl = 0, double ch = 180)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   bool on[GNL] = {true, true, true, true};
   const char *sn[3] = {"generated truth", "accepted truth", "reconstructed"};
   auto *c = new TCanvas("cKinPng", "kin", 1400, 950);
   c->SetRightMargin(0.13);
   c->SetLeftMargin(0.10);
   c->SetLogz();
   auto *h = new TH2D("hKinPng",
                      Form("^{10}Be(t,p)^{12}Be   %s   %s;#theta_{lab} [deg];KE_{p} [MeV]",
                           samp == 0 ? "(field-independent)" : cfg.Data(), sn[samp]),
                      360, 0, 180, 275, 0, 55);
   h->SetDirectory(nullptr);
   kin_fill(root, outDir, cfg, samp, on, cl, ch, h, nullptr);
   h->Draw("colz");
   auto *leg = new TLegend(0.55, 0.66, 0.86, 0.90);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.026);
   for (int l = 0; l < GNL; ++l) {
      auto *g = new TGraph();
      int n = 0;
      for (double a = 0.5; a <= 179.5; a += 0.25) {
         double th, ke;
         if (!g_fwd(GM4 + GLEX[l], a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
         g->SetPoint(n++, th * TMath::RadToDeg(), ke);
      }
      g->SetLineColor(GLCOL[l]);
      g->SetLineWidth(3);
      g->Draw("L SAME");
      leg->AddEntry(g, GLJP[l], "l");
   }
   leg->Draw();
   auto *gs = new TGraph();
   auto *tx = new TLatex();
   tx->SetTextSize(0.020);
   tx->SetTextColor(kGray + 3);
   int n = 0;
   for (double a : {5., 10., 20., 30., 45., 60., 90., 120., 150.}) {
      double th, ke;
      if (!g_fwd(GM4, a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
      gs->SetPoint(n++, th * TMath::RadToDeg(), ke);
      tx->DrawLatex(th * TMath::RadToDeg() + 1.5, ke + 1.0, Form("%.0f#circ", a));
   }
   gs->SetMarkerStyle(29);
   gs->SetMarkerSize(1.5);
   gs->SetMarkerColor(kGray + 3);
   gs->Draw("P SAME");
   auto *lb = new TLatex();
   lb->SetNDC();
   lb->SetTextSize(0.022);
   lb->SetTextColor(kGray + 3);
   lb->DrawLatex(0.14, 0.16, "stars: #theta_{cm} (DWBA) -- runs BACKWARDS across this plot");
   TString f = TString::Format("%s/kin_%s_s%d.png", outDir.Data(), cfg.Data(), samp);
   c->SaveAs(f);
   printf("wrote %s  (%.0f entries)\n", f.Data(), h->GetEntries());
}
