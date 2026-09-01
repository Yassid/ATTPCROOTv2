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
   TGCheckButton *fLev[GNL]{}, *fCurves{nullptr}, *fCmLines{nullptr}, *fLogZ{nullptr}, *fVtxE{nullptr};
   TGNumberEntry *fCmLo{nullptr}, *fCmHi{nullptr};
   TGNumberEntry *fNx{nullptr}, *fXlo{nullptr}, *fXhi{nullptr};
   TGNumberEntry *fNy{nullptr}, *fYlo{nullptr}, *fYhi{nullptr};
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
   void FitRange();
   void SetPeakWindow();
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

   // BINNING AND RANGE. Both, together: rebinning without being able to move the range is half a
   // control -- the interesting structure here (the four loci separating at forward angle, the
   // 2.109/2.251 pair merging at backward angle) only shows when the window follows the bins in.
   ctl->AddFrame(new TGLabel(ctl, "binning and range"), new TGLayoutHints(kLHintsLeft, 4, 4, 10, 2));
   auto mkrow = [&](const char *lab, TGNumberEntry *&n_, TGNumberEntry *&lo_, TGNumberEntry *&hi_, double nd,
                    double lod, double hid, int idbase) {
      auto *row = new TGHorizontalFrame(ctl);
      row->AddFrame(new TGLabel(row, lab), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 2, 0, 0));
      n_ = new TGNumberEntry(row, nd, 4, idbase, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
      lo_ = new TGNumberEntry(row, lod, 4, idbase + 1, TGNumberFormat::kNESRealOne);
      hi_ = new TGNumberEntry(row, hid, 4, idbase + 2, TGNumberFormat::kNESRealOne);
      n_->GetNumberEntry()->SetToolTipText("number of bins");
      lo_->GetNumberEntry()->SetToolTipText("axis minimum");
      hi_->GetNumberEntry()->SetToolTipText("axis maximum");
      row->AddFrame(n_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
      row->AddFrame(lo_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
      row->AddFrame(hi_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
      ctl->AddFrame(row, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));
   };
   ctl->AddFrame(new TGLabel(ctl, "        bins     min      max"), new TGLayoutHints(kLHintsLeft, 4, 4, 0, 0));
   mkrow("#lab", fNx, fXlo, fXhi, 360, 0, 180, 500);
   mkrow("KE ", fNy, fYlo, fYhi, 275, 0, 55, 510);
   auto *bPeak = new TGTextButton(ctl, "set #theta_{cm} = transfer peak (2-45)", 405);
   ctl->AddFrame(bPeak, new TGLayoutHints(kLHintsExpandX, 8, 8, 6, 2));
   auto *bFit = new TGTextButton(ctl, "fit range to data", 404);
   ctl->AddFrame(bFit, new TGLayoutHints(kLHintsExpandX, 8, 8, 4, 6));

   fCurves = new TGCheckButton(ctl, "analytic level curves", 300);
   fCurves->SetOn(kTRUE);
   ctl->AddFrame(fCurves, new TGLayoutHints(kLHintsLeft, 12, 4, 6, 1));
   fCmLines = new TGCheckButton(ctl, "iso-#theta_{cm} markers", 301);
   fCmLines->SetOn(kTRUE);
   ctl->AddFrame(fCmLines, new TGLayoutHints(kLHintsLeft, 12, 4, 1, 1));
   fLogZ = new TGCheckButton(ctl, "log z", 302);
   fLogZ->SetOn(kTRUE);
   ctl->AddFrame(fLogZ, new TGLayoutHints(kLHintsLeft, 12, 4, 1, 1));
   // Worth having as a TOGGLE rather than just switching it on: flipping it is the clearest
   // demonstration in the whole study of what the vertex beam-energy correction is worth
   // (E_x core width 0.230 -> 0.093 MeV on the ground state at theta_cm 2-45).
   fVtxE = new TGCheckButton(ctl, "vertex E_{beam} correction", 303);
   fVtxE->SetOn(kTRUE);
   ctl->AddFrame(fVtxE, new TGLayoutHints(kLHintsLeft, 12, 4, 1, 8));

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
   fVtxE->Connect("Clicked()", "KinGui", this, "Redraw()");
   bFit->Connect("Clicked()", "KinGui", this, "FitRange()");
   bPeak->Connect("Clicked()", "KinGui", this, "SetPeakWindow()");
   for (TGNumberEntry *e : {fNx, fXlo, fXhi, fNy, fYlo, fYhi, fCmLo, fCmHi})
      e->Connect("ValueSet(Long_t)", "KinGui", this, "Redraw()");
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

/// Ex from (theta_lab, KE) at a GIVEN beam energy, for the ground-state residual mass.
static double g_exAt(double Eb, double th, double Ke)
{
   double Et1 = Eb + GM1, Et3 = Ke + GM3;
   double s = GM1 * GM1 + GM2 * GM2 + 2 * GM2 * Et1;
   double uu = GM2 * GM2 + GM3 * GM3 - 2 * GM2 * Et3;
   double a = (std::cos(th) * g_om2(s, GM1 * GM1, GM2 * GM2) * g_om2(uu, GM2 * GM2, GM3 * GM3) -
               (s - GM1 * GM1 - GM2 * GM2) * (GM2 * GM2 + GM3 * GM3 - uu)) / (2 * GM2 * GM2) + s + uu - GM2 * GM2;
   return a > 0 ? std::sqrt(a) - GM4 : NAN;
}

/// THE BEAM ENERGY AT THE VERTEX, and why this function exists.
///
/// The beam loses 5.3 MeV crossing the metre of T2, so the energy at the reaction is a function of
/// where the reaction happened. Holding it constant costs dEx/dE_beam x (the spread of vertex
/// energies) = 0.13 x 1.4 MeV of excitation-energy resolution, which is most of what this channel
/// has. Measured on the ground state at theta_cm 2-45: the E_x core width is 0.230 MeV with a
/// constant beam energy and 0.093 MeV with the vertex one.
///
/// The GUI's E_x projection originally plotted the exres tree's `exReco` branch, which
/// ex_res_C14_hf.C fills with ONE constant Ebeam for every event -- so the projection looked far
/// worse than tp_spectrum_Be10.C's for a reason that had nothing to do with the detector. This is
/// the same profile the spectrum macro fits: E_beam solved from truth per event, then a quadratic
/// in the vertex z.
static void g_ebeamProfile(TTree *t, double lvlEx, const char *bTh, const char *bKe, const char *bZ, TF1 &fEb)
{
   double th, ke, z;
   t->SetBranchAddress(bTh, &th);
   t->SetBranchAddress(bKe, &ke);
   t->SetBranchAddress(bZ, &z);
   std::vector<double> eb, zz;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      double lo = 80., hi = 130.;
      double thr = th * TMath::DegToRad();
      double flo = g_exAt(lo, thr, ke) - lvlEx, fhi = g_exAt(hi, thr, ke) - lvlEx;
      if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
      for (int it = 0; it < 60; ++it) {
         double m = 0.5 * (lo + hi), fm = g_exAt(m, thr, ke) - lvlEx;
         if (std::isnan(fm)) break;
         if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
      }
      double e = 0.5 * (lo + hi);
      if (e > 85 && e < 125) { eb.push_back(e); zz.push_back(z); }
   }
   fEb.SetParameters(GEB, -0.006, 0.);
   if (eb.size() >= 100) {
      TGraph g((int)eb.size(), zz.data(), eb.data());
      g.Fit(&fEb, "QN");
   }
}

/// Fill the 2D (and, if asked, the E_x spectrum) for one selection. A FREE function, not a member,
/// so the batch renderer below draws from exactly the same code the GUI does -- a figure that came
/// from a second copy of this loop would be a figure of something slightly different.
///   samp 0 = generated truth (from the cache), 1 = accepted truth, 2 = reconstructed
static void kin_fill(const TString &root, const TString &outDir, const TString &cfg, int samp, const bool *levOn,
                     double cl, double ch, TH2D *h, TH1D *hex, bool vtxE = true)
{
   for (int l = 0; l < GNL; ++l) {
      if (!levOn[l]) continue;
      if (samp == 0) {
         TFile *f = TFile::Open(outDir + "/kin_truth_Be10.root");
         if (!f || f->IsZombie()) continue;
         TTree *t = (TTree *)f->Get(Form("truth_%s", GLVN[l]));
         if (!t) { f->Close(); continue; }
         TF1 fEb("fEbT", "[0]+[1]*x+[2]*x*x", 0, 1000);
         if (vtxE) g_ebeamProfile(t, GLEX[l], "th", "ke", "z", fEb);
         double th, ke, z;
         t->SetBranchAddress("th", &th);
         t->SetBranchAddress("ke", &ke);
         t->SetBranchAddress("z", &z);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double ex, cm;
            g_inv(th * TMath::DegToRad(), ke, ex, cm);
            if (!(cm >= cl && cm < ch)) continue;
            if (h) h->Fill(th, ke);
            if (hex) hex->Fill(vtxE ? g_exAt(fEb.Eval(z), th * TMath::DegToRad(), ke) : ex);
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
         TF1 fEb("fEbR", "[0]+[1]*x+[2]*x*x", 0, 1000);
         if (vtxE) g_ebeamProfile(t, GLEX[l], "thTrue", "keTrue", "zTrue", fEb);
         double thT, keT, thR, keR, cmT, exR, zT, zR;
         t->SetBranchAddress("thTrue", &thT);
         t->SetBranchAddress("keTrue", &keT);
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keReco", &keR);
         t->SetBranchAddress("cmTrue", &cmT);
         t->SetBranchAddress("exReco", &exR);
         t->SetBranchAddress("zTrue", &zT);
         t->SetBranchAddress("zReco", &zR);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!(cmT >= cl && cmT < ch)) continue;
            if (samp == 1) {
               // ACCEPTED TRUTH. This used to Fill(GLEX[l]) -- the exact level energy, i.e. a
               // delta spike, which is not a spectrum of anything. What it should show is the
               // RESOLUTION FLOOR: the truth angle and energy pushed through the same inversion
               // the analysis uses, so the only width left is what the beam-energy treatment and
               // the vertex spread put there. That is the floor no fitter can beat.
               if (h) h->Fill(thT, keT);
               if (hex) hex->Fill(vtxE ? g_exAt(fEb.Eval(zT), thT * TMath::DegToRad(), keT)
                                       : g_exAt(GEB, thT * TMath::DegToRad(), keT));
            } else {
               if (h) h->Fill(thR, keR);
               // exReco is the CONSTANT-Ebeam value the producer stored; recompute at the
               // reconstructed vertex when the correction is on.
               if (hex) hex->Fill(vtxE ? g_exAt(fEb.Eval(zR), thR * TMath::DegToRad(), keR) : exR);
            }
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
            fCmHi->GetNumber(), h, hex, fVtxE->IsOn());
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

   int nx = (int)fNx->GetNumber(), ny = (int)fNy->GetNumber();
   double xlo = fXlo->GetNumber(), xhi = fXhi->GetNumber();
   double ylo = fYlo->GetNumber(), yhi = fYhi->GetNumber();
   // A reversed or empty axis makes TH2D throw and takes the GUI down with it, so clamp rather
   // than trust the entry boxes -- a user halfway through typing "180" has momentarily asked for
   // a max of 1.
   if (!(xhi > xlo)) xhi = xlo + 1;
   if (!(yhi > ylo)) yhi = ylo + 1;
   if (nx < 1) nx = 1;
   if (ny < 1) ny = 1;
   fH = new TH2D("hKin", Form("^{10}Be(t,p)^{12}Be   %s   %s;#theta_{lab} [deg];KE_{p} [MeV]",
                              samp == 0 ? "(field-independent)" : cfg.Data(), sn[samp]),
                 nx, xlo, xhi, ny, ylo, yhi);
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
   fBar->SetText(Form("%s | %s | %.0f entries | bins %dx%d (%.2f deg x %.3f MeV)", cfg.Data(), sn[samp],
                      fH->GetEntries(), nx, ny, (xhi - xlo) / nx, (yhi - ylo) / ny),
                 0);
}

void KinGui::OnMove(Int_t, Int_t px, Int_t py, TObject *)
{
   TCanvas *c = fEc->GetCanvas();
   double x = c->AbsPixeltoX(px), y = c->AbsPixeltoY(py);
   double th = c->PadtoX(x), ke = c->PadtoY(y);
   if (!fH) return;
   if (th < fH->GetXaxis()->GetXmin() || th > fH->GetXaxis()->GetXmax() || ke < fH->GetYaxis()->GetXmin() ||
       ke > fH->GetYaxis()->GetXmax())
      return;
   double ex, cm;
   g_inv(th * TMath::DegToRad(), ke, ex, cm);
   fBar->SetText(Form("theta_lab %.2f deg   KE_p %.3f MeV", th, ke), 1);
   if (std::isnan(ex)) fBar->SetText("outside the kinematic locus", 2);
   else fBar->SetText(Form("=>  E_x %+.3f MeV   theta_cm %.1f deg", ex, std::isnan(cm) ? 0.0 : cm), 2);
}

void KinGui::ProjectEx()
{
   auto *c = new TCanvas(Form("cEx%d", (int)gRandom->Integer(100000)), "Ex projection", 1000, 650);
   auto *h = new TH1D("hExProj", Form("E_{x} projection -- %s beam energy;E_{x}(^{12}Be) [MeV];counts",
                                      fVtxE->IsOn() ? "VERTEX" : "constant"),
                      300, -1.5, 4.5);
   Fill(nullptr, h);
   h->SetLineColor(kBlack);
   h->Draw();
   // Quote the width on the plot. Without it the two beam-energy treatments look like two
   // pictures rather than a factor of 2.5.
   {
      double q[3] = {0.25, 0.5, 0.75}, v[3];
      h->GetQuantiles(3, v, q);
      double n = h->GetEntries(), tail = 0;
      for (int b = 1; b <= h->GetNbinsX(); ++b)
         if (std::fabs(h->GetBinCenter(b) - v[1]) > 1.0) tail += h->GetBinContent(b);
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.033);
      tx->DrawLatex(0.14, 0.87, Form("IQR/1.349 = %.3f MeV   (%s E_{beam})", (v[2] - v[0]) / 1.349,
                                     fVtxE->IsOn() ? "vertex" : "constant"));
      // The TAIL FRACTION is the number that explains a bad-looking projection, and it is almost
      // always the angular window rather than anything else: 3 % beyond +-1 MeV over theta_cm
      // 2-45, but 25 % over 2-180.
      tx->DrawLatex(0.14, 0.82, Form("beyond #pm1 MeV: %.1f %%   (#theta_{cm} %.0f-%.0f#circ)",
                                     n > 0 ? 100 * tail / n : 0.0, fCmLo->GetNumber(), fCmHi->GetNumber()));
      if (fCmHi->GetNumber() > 50) {
         tx->SetTextColor(kRed + 1);
         tx->DrawLatex(0.14, 0.77, "window includes #theta_{cm} > 50#circ: fast forward protons,");
         tx->DrawLatex(0.14, 0.73, "#sigma(KE)/KE 3-6 %, irreducible. Use the transfer peak.");
      }
   }
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

/// THE ANGULAR WINDOW IS THE OTHER HALF OF THE RESOLUTION, and it is easy to forget because the
/// default (everything) looks like the honest choice. Measured on the ground state, E_x width with
/// the vertex beam energy already applied:
///
///     theta_cm     const Eb   vertex Eb   floor    median KE   sigma(KE)/KE
///     2-20           0.195      0.077     0.020      4.2 MeV       0.73 %
///     20-45          0.241      0.097     0.021      7.3 MeV       0.67 %
///     45-90          0.653      0.601     0.026     18.4 MeV       3.42 %
///     90-180         1.364      1.369     0.032     37.4 MeV       5.66 %
///
/// The vertex correction is worth a factor 2.5 backward and NOTHING forward -- at theta_cm 90-180
/// it is fractionally worse. The reason is in the last column: a fast forward proton has a helix
/// radius comparable to the chamber and traverses a short arc, so its curvature, hence its
/// momentum, is badly determined. 5.7 % of 37 MeV is 2.1 MeV, and dEx/dKE ~ 0.7 turns that into
/// 1.4 MeV of excitation energy. No correction reaches it: the floor is 0.03 MeV, so it is all
/// detector. Opening the window to 0-180 puts 43 % of the sample in that regime, which is exactly
/// the tail that makes this projection look worse than tp_spectrum_Be10.C's.
void KinGui::SetPeakWindow()
{
   fCmLo->SetNumber(2);
   fCmHi->SetNumber(45);
   Redraw();
}

/// Snap the axes to the data actually selected, keeping the bin WIDTH the user has chosen rather
/// than the bin COUNT -- zooming in should give more bins over the same structure, not the same
/// number of wider ones stretched over less range.
void KinGui::FitRange()
{
   if (!fH) return;
   double wx = (fH->GetXaxis()->GetXmax() - fH->GetXaxis()->GetXmin()) / fH->GetNbinsX();
   double wy = (fH->GetYaxis()->GetXmax() - fH->GetYaxis()->GetXmin()) / fH->GetNbinsY();
   int x1 = fH->GetXaxis()->GetFirst(), x2 = fH->GetXaxis()->GetLast();
   int y1 = fH->GetYaxis()->GetFirst(), y2 = fH->GetYaxis()->GetLast();
   double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9;
   for (int i = x1; i <= x2; ++i)
      for (int j = y1; j <= y2; ++j)
         if (fH->GetBinContent(i, j) > 0) {
            xlo = std::min(xlo, fH->GetXaxis()->GetBinLowEdge(i));
            xhi = std::max(xhi, fH->GetXaxis()->GetBinUpEdge(i));
            ylo = std::min(ylo, fH->GetYaxis()->GetBinLowEdge(j));
            yhi = std::max(yhi, fH->GetYaxis()->GetBinUpEdge(j));
         }
   if (!(xhi > xlo) || !(yhi > ylo)) { fBar->SetText("nothing to fit the range to", 0); return; }
   double px = 0.02 * (xhi - xlo), py = 0.02 * (yhi - ylo);
   xlo -= px; xhi += px; ylo -= py; yhi += py;
   fXlo->SetNumber(xlo); fXhi->SetNumber(xhi);
   fYlo->SetNumber(ylo); fYhi->SetNumber(yhi);
   fNx->SetNumber(std::max(1.0, std::floor((xhi - xlo) / wx)));
   fNy->SetNumber(std::max(1.0, std::floor((yhi - ylo) / wy)));
   Redraw();
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
                  TString outDir = "plots", double cl = 0, double ch = 180, int nx = 360, double xlo = 0,
                  double xhi = 180, int ny = 275, double ylo = 0, double yhi = 55)
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
                      nx, xlo, xhi, ny, ylo, yhi);
   h->SetDirectory(nullptr);
   kin_fill(root, outDir, cfg, samp, on, cl, ch, h, nullptr, true);
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
