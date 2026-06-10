/// @file replot_pd_a1975.C
/// @brief Regenerate the 16C(p,d)15C 4-panel from the cached ntuple (fast, no file
///        re-read), with the g.s.->0 energy calibration applied.
///
/// Same 4 panels as ex_pd_a1975.C — 15C Ex (calibrated), Ex vs theta_cm, Brho vs
/// theta_lab, PID plane — but rebuilt from deuteron_kin.root in seconds. Ex is
/// shifted by exShift (default -0.615; the cache stores raw Ex).
///
///   root -b -q 'pd/replot_pd_a1975.C'                 // calibrated
///   root -b -q 'pd/replot_pd_a1975.C(0)'              // raw

void replot_pd_a1975(double exShift = -0.615, TString cut = "", TString cacheFile = "deuteron_kin.root",
                     TString outPng = "pd/plots/ex_pd_spectrum.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   if (!t) {
      printf("ERROR: ntuple 'dk' not in %s\n", cacheFile.Data());
      return;
   }

   TString exC = Form("ex+(%g)", exShift); // calibrated Ex expression
   const char *tit = (exShift != 0) ? "calibrated" : "raw";

   TH1F *hex = new TH1F("hex", Form("^{15}C excitation energy (p,d), %s;E_{x}(^{15}C) [MeV];deuteron candidates", tit),
                        200, -10, 20);
   TH2F *hexcm =
      new TH2F("hexcm", "E_{x}(^{15}C) vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -10, 20);
   TH2F *hbt = new TH2F("hbt", "B#rho vs #theta_{lab} (deuteron cand.);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200,
                        0, 2.0);
   TH2F *hpid = new TH2F("hpid", "PID plane (deuteron cand.);#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 2.0);

   t->Project("hex", exC, cut);
   t->Project("hexcm", exC + ":thcm", cut);
   t->Project("hbt", "brho:theta", cut);
   t->Project("hpid", "brho:sqrtdedx", cut);
   printf("replotted %.0f candidates (exShift=%.3f%s)\n", hex->GetEntries(), exShift,
          cut.Length() ? Form(", cut: %s", cut.Data()) : "");

   TCanvas *c = new TCanvas("c", "ex_pd", 1550, 1040);
   c->Divide(2, 2);
   c->cd(1);
   hex->SetLineColor(kBlue + 1);
   hex->SetLineWidth(2);
   hex->Draw("hist");
   // known 15C levels
   double lev[] = {0.0, 0.740, 3.103, 4.220, 4.657, 5.866};
   const char *lbl[] = {"g.s.", "0.74", "3.10", "4.22", "4.66", "5.87"};
   hex->GetXaxis()->SetRangeUser(-3, 11);
   double ymax = hex->GetMaximum();
   for (int i = 0; i < 6; ++i) {
      TLine *l = new TLine(lev[i], 0, lev[i], ymax);
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
      TLatex *tx = new TLatex(lev[i] + 0.05, ymax * (0.97 - 0.06 * (i % 2)), lbl[i]);
      tx->SetTextColor(kRed + 1);
      tx->SetTextSize(0.03);
      tx->Draw();
   }
   c->cd(2);
   hexcm->Draw("colz");
   c->cd(3);
   hbt->Draw("colz");
   c->cd(4);
   hpid->Draw("colz");
   c->SaveAs(outPng);
   printf("saved %s\n", outPng.Data());
}
