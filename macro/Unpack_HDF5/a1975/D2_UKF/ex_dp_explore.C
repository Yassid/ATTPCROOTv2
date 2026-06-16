/// @file ex_dp_explore.C
/// @brief INTERACTIVE 17C(d,p) Ex explorer. Loads the cached proton kinematics once
///        and lets you re-plot any proton theta_lab window live at the ROOT prompt —
///        the canvas stays interactive (zoom, hover, save) and `exr(...)` redraws it
///        instantly. No re-fit / no re-read of the genfitter files.
///
/// Run it WITHOUT -b and WITHOUT -q so ROOT keeps the canvas and the prompt:
///
///   cd .../a1975/D2_UKF
///   source .../build/config.sh
///   root -l ex_dp_explore.C                 // opens canvas at the default backward window
///
/// Then at the `root [n]` prompt, redraw any window/cut live:
///   exr(90,180)        // backward (d,p) signal
///   exr(0,90)          // forward
///   exr(90,140, 5)     // + chi2/ndf < 5
///   exr(90,180, 0, -7.5)   // slide Ex (g.s.->0) once Ebeam is calibrated
///   .q                 // quit
///
/// Args of exr(): thLoDeg, thHiDeg, chi2Max (0=off), exShift.
/// A 2D TCanvas works under WSLg (it is ROOT's 3D/Eve viewers that crash there).

namespace {
TNtuple *gPK = nullptr; // the cached ntuple, loaded once
TCanvas *gExC = nullptr; // the persistent interactive canvas
TString gCache;
} // namespace

/// Redraw the 4-panel view for a theta_lab window. Callable repeatedly at the prompt.
void exr(double thLoDeg = 90.0, double thHiDeg = 180.0, double chi2Max = 0.0, double exShift = 0.0)
{
   if (!gPK) {
      printf("ntuple not loaded — run ex_dp_explore.C first (root -l ex_dp_explore.C)\n");
      return;
   }
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString cut = Form("theta>=%g&&theta<=%g", thLoDeg, thHiDeg);
   if (chi2Max > 0)
      cut += Form("&&chi2ndf<%g", chi2Max);
   TString exExpr = (exShift != 0.0) ? Form("(ex+%g)", exShift) : TString("ex");
   long nb = (long)gPK->GetEntries(cut);
   TString tag = (thLoDeg >= 90) ? "BACKWARD" : (thHiDeg <= 90 ? "FORWARD" : "");

   // drop any previous histos so re-draws don't accumulate / warn
   for (auto nm : {"hexX", "hexcmX", "hbtX", "hpidX"})
      if (auto *o = gDirectory->Get(nm))
         delete o;

   auto *hex = new TH1F(
      "hexX",
      Form("^{17}C E_{x} (d,p) %s #theta_{lab} %.0f-%.0f#circ  [%ld cand%s];E_{x}(^{17}C) [MeV];proton candidates",
           tag.Data(), thLoDeg, thHiDeg, nb, chi2Max > 0 ? Form(", chi2<%g", chi2Max) : ""),
      200, -10, 20);
   auto *hexcm = new TH2F("hexcmX", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -10, 20);
   auto *hbt = new TH2F("hbtX", "B#rho vs #theta_{lab};#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200, 0, 1.2);
   auto *hpid = new TH2F("hpidX", "PID plane;#sqrt{dEdx};B#rho [T m]", 200, 0, 50, 200, 0, 1.2);

   gPK->Draw(exExpr + ">>hexX", cut, "goff");
   gPK->Draw(exExpr + ":thcm>>hexcmX", cut, "goff");
   gPK->Draw("brho:theta>>hbtX", cut, "goff");
   gPK->Draw("brho:sqrtdedx>>hpidX", cut, "goff");

   if (!gExC || !gROOT->GetListOfCanvases()->FindObject(gExC))
      gExC = new TCanvas("exExplore", "17C(d,p) Ex explorer", 1550, 1040);
   gExC->Clear();
   gExC->Divide(2, 2);
   gExC->cd(1);
   hex->Draw();
   auto *l0 = new TLine(0, 0, 0, hex->GetMaximum());
   l0->SetLineColor(kRed);
   l0->SetLineStyle(2);
   l0->Draw();
   gExC->cd(2);
   hexcm->Draw("colz");
   gExC->cd(3);
   hbt->Draw("colz");
   gExC->cd(4);
   hpid->Draw("colz");
   gExC->Update();

   printf("theta_lab %.0f-%.0f: %ld candidates  Ex peak=%.2f MeV  mean=%.2f\n", thLoDeg, thHiDeg, nb,
          hex->GetBinCenter(hex->GetMaximumBin()), hex->GetMean());
}

/// Entry point: load the cache once, draw an initial window, print the usage hint.
void ex_dp_explore(double thLoDeg = 90.0, double thHiDeg = 180.0, const char *cache = "proton_kin_dp.root")
{
   if (!gPK || gCache != cache) {
      TFile *f = TFile::Open(cache);
      gPK = f ? (TNtuple *)f->Get("pk") : nullptr;
      gCache = cache;
   }
   if (!gPK) {
      printf("no ntuple 'pk' in %s — run ex_dp_a1975.C first\n", cache);
      return;
   }
   exr(thLoDeg, thHiDeg);
   printf("\n\033[1;33m--- interactive: redraw any window at the prompt ---\033[0m\n");
   printf("  exr(90,180)        backward (d,p) signal\n");
   printf("  exr(0,90)          forward\n");
   printf("  exr(90,140, 5)     + chi2/ndf<5\n");
   printf("  exr(90,180, 0, -7.5)   slide Ex (g.s.->0)\n");
   printf("  .q                 quit\n");
}
