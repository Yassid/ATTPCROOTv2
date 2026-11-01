/// @file plot_pid_spyral_vs_best_dp.C
/// @brief Spyral with run + pad gain (sqp_g vs br) against the adopted PID (3D s1_g vs bA below bSwitch, br above) on
///        the SAME tracks (2026-11-01). Top: the two planes, Brho on y, same frame. Bottom: ln-mass per Brho slice
///        (each method sliced by its own Brho, each curve shifted so ITS deuteron mode in that slice is at 0), log y.
///   root -l -b -q 'plot_pid_spyral_vs_best_dp.C("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root","plots/pid_spyral_vs_best.png")'
void plot_pid_spyral_vs_best_dp(TString in, TString png, double bSwitch = 0.45,
                                TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&br>0&&sqp_g>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   const char *nm[2] = {"Spyral, run + pad gain", "adopted: 3D dE/dx, hit-refit / Spyral B#rho"};
   TString sq[2] = {"sqp_g", "s1_g"}, B[2] = {"br", Form("(bA<%g?bA:br)", bSwitch)};
   double med[2], yhi[2], dref[2];
   for (int m = 0; m < 2; ++m) {
      w->Draw(sq[m], sel, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      med[m] = x[x.size() / 2]; yhi[m] = x[x.size() * 995 / 1000] * 1.1;
      // reference deuteron position: mode at Brho 0.15-0.35 (deuteron dominated)
      TH1D hr(Form("hr%d", m), "", 600, -4, 2);
      w->Draw(Form("log(%s/%g*%s)>>hr%d", sq[m].Data(), med[m], B[m].Data(), m),
              sel + Form("&&%s>0.15&&%s<0.35", B[m].Data(), B[m].Data()), "goff");
      hr.Smooth(2);
      dref[m] = hr.GetBinCenter(hr.GetMaximumBin());
   }
   auto *c = new TCanvas("c", "c", 2000, 1400);
   auto *top = new TPad("top", "", 0, 0.45, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.45); bot->Draw();
   top->Divide(2, 1);
   for (int m = 0; m < 2; ++m) {
      top->cd(m + 1);
      gPad->SetLogz(); gPad->SetRightMargin(0.12);
      // dE/dx in units of the median so the two frames match
      auto *h2 = new TH2D(Form("p%d", m), Form("%s;#sqrt{dE/dx} / median;B#rho [T m]", nm[m]), 400, 0, 2.6, 400, 0, 1.2);
      w->Draw(Form("%s:%s/%g>>p%d", B[m].Data(), sq[m].Data(), med[m], m), sel, "goff");
      h2->Draw("colz");
   }
   bot->Divide(5, 1);
   const double bl[] = {0.10, 0.20, 0.30, 0.45, 0.60, 0.80};
   const int col[2] = {kRed + 1, kBlue + 2};
   for (int b = 0; b < 5; ++b) {
      bot->cd(b + 1);
      gPad->SetLogy(); gPad->SetGridx();
      TH1D *h[2];
      double ymax = 0;
      for (int m = 0; m < 2; ++m) {
         TString L = Form("log(%s/%g*%s)", sq[m].Data(), med[m], B[m].Data());
         TString cut = sel + Form("&&%s>%g&&%s<%g", B[m].Data(), bl[b], B[m].Data(), bl[b + 1]);
         TH1D tmp(Form("t%d%d", b, m), "", 600, -4, 2);
         w->Draw(L + Form(">>t%d%d", b, m), cut, "goff");
         // deuteron mode: highest smoothed bin within -0.25..+0.3 of the reference d position (global max can be p/4He)
         TH1D s(tmp); s.Smooth(2);
         int i0 = s.FindBin(dref[m] - 0.25), i1 = s.FindBin(dref[m] + 0.3), bd = i0;
         for (int i = i0; i <= i1; ++i) if (s.GetBinContent(i) > s.GetBinContent(bd)) bd = i;
         double d = s.GetBinCenter(bd);
         h[m] = new TH1D(Form("h%d%d", b, m), Form("B#rho %.2f-%.2f T m;ln mass - ln m_{d};tracks", bl[b], bl[b + 1]), 200, -1.2, 1.2);
         w->Draw(Form("%s-(%g)>>h%d%d", L.Data(), d, b, m), cut, "goff");
         h[m]->SetLineColor(col[m]); h[m]->SetLineWidth(2);
         ymax = std::max(ymax, h[m]->GetMaximum());
      }
      h[0]->SetMaximum(3 * ymax); h[0]->SetMinimum(0.8);
      h[0]->Draw("hist"); h[1]->Draw("hist same");
      const double lm[] = {std::log(938.272 / 1875.613), 0, std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)};
      for (double v : lm) { auto *l = new TLine(v, 0.8, v, 3 * ymax); l->SetLineStyle(3); l->Draw(); }
      if (b == 0) {
         auto *lg = new TLegend(0.12, 0.78, 0.7, 0.89);
         lg->AddEntry(h[0], "Spyral run+pad", "l"); lg->AddEntry(h[1], "adopted", "l"); lg->Draw();
      }
   }
   c->SaveAs(png);
}
