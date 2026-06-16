/// @file ex_dp_range.C
/// @brief 17C(d,p) Ex spectrum in a selectable PROTON theta_lab range, from the cached
///        proton kinematics ntuple (proton_kin_dp.root, produced by ex_dp_a1975.C).
///        Same 4-panel view as ex_dp_a1975.C but with a theta_lab window so you can
///        isolate, e.g., the BACKWARD (d,p) signal (theta_lab 90-180) from the
///        forward beam-like junk. No re-fit / no re-read of the genfitter files.
///
///   root -l -b -q 'ex_dp_range.C(90,180)'                 // backward only (the (d,p) signal)
///   root -l -b -q 'ex_dp_range.C(0,90)'                   // forward only
///   root -l -b -q 'ex_dp_range.C(90,180,5.0,"proton_kin_dp.root",-7.5)'  // + chi2<5 + g.s. shift
///
/// Args: thLoDeg, thHiDeg (proton theta_lab window, deg); chi2Max (extra chi2/ndf cut,
/// <=0 = none); cache (ntuple file); exShift (slide Ex, e.g. to put g.s. at 0).

void ex_dp_range(double thLoDeg = 90.0, double thHiDeg = 180.0, double chi2Max = 0.0,
                 const char *cache = "proton_kin_dp.root", double exShift = 0.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *f = TFile::Open(cache);
   TNtuple *t = f ? (TNtuple *)f->Get("pk") : nullptr;
   if (!t) {
      printf("no ntuple 'pk' in %s — run ex_dp_a1975.C first\n", cache);
      return;
   }

   // cut string: theta_lab window (+ optional chi2). exShift applied to ex via (ex+shift).
   TString cut = Form("theta>=%g&&theta<=%g", thLoDeg, thHiDeg);
   if (chi2Max > 0)
      cut += Form("&&chi2ndf<%g", chi2Max);
   TString exExpr = (exShift != 0.0) ? Form("(ex+%g)", exShift) : TString("ex");

   long nb = (long)t->GetEntries(cut);
   TString tag = (thLoDeg >= 90) ? "BACKWARD" : (thHiDeg <= 90 ? "FORWARD" : "");

   TH1F *hex = new TH1F(
      "hexr",
      Form("^{17}C E_{x} (d,p), %s #theta_{lab} %.0f-%.0f#circ  [%ld cand];E_{x}(^{17}C) [MeV];proton candidates",
           tag.Data(), thLoDeg, thHiDeg, nb),
      200, -10, 20);
   TH2F *hexcm = new TH2F("hexcmr", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -10, 20);
   TH2F *hbt = new TH2F("hbtr", "B#rho vs #theta_{lab};#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200, 0, 1.2);
   TH2F *hpid = new TH2F("hpidr", "PID plane;#sqrt{dEdx};B#rho [T m]", 200, 0, 50, 200, 0, 1.2);

   t->Draw(exExpr + ">>hexr", cut, "goff");
   t->Draw(exExpr + ":thcm>>hexcmr", cut, "goff");
   t->Draw("brho:theta>>hbtr", cut, "goff");
   t->Draw("brho:sqrtdedx>>hpidr", cut, "goff");

   TCanvas *c = new TCanvas("cr", "ex_dp_range", 1550, 1040);
   c->Divide(2, 2);
   c->cd(1);
   hex->Draw();
   TLine *l0 = new TLine(0, 0, 0, hex->GetMaximum());
   l0->SetLineColor(kRed);
   l0->SetLineStyle(2);
   l0->Draw();
   c->cd(2);
   hexcm->Draw("colz");
   c->cd(3);
   hbt->Draw("colz");
   c->cd(4);
   hpid->Draw("colz");

   TString png = Form("plots/ex_dp_range_%.0f_%.0f.png", thLoDeg, thHiDeg);
   c->SaveAs(png);
   printf("theta_lab %.0f-%.0f: %ld candidates  Ex peak=%.2f MeV  mean=%.2f  -> %s\n", thLoDeg, thHiDeg, nb,
          hex->GetBinCenter(hex->GetMaximumBin()), hex->GetMean(), png.Data());
}
