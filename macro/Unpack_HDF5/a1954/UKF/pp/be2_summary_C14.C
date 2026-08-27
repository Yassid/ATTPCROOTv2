/// @file be2_summary_C14.C
/// @brief The 7.012 B(E2) under each candidate potential, against the electromagnetic value.
///
/// Each potential is used CONSISTENTLY -- its own luminosity from the elastic AND its own DWBA --
/// so the partial cancellation between the two is included (see exc_omp_consistent_C14.C). The
/// error bar on each is the 6.903 blend band (f903 = 0 to 0.5, which the data cannot resolve at
/// 0.76 sigma separation) combined with the +-13% potential spread.
///
/// The literature value is ENSDF's B(E2; 2+ -> 0+) = 1.8(3) W.u. from the 9.0 fs lifetime,
/// converted UP by the factor (2Jf+1)/(2Ji+1) = 5, because Ptolemy's BELX is the excitation
/// strength while tabulated lifetimes give the decay one. That factor alone looked like most of
/// an order of magnitude before it was spotted.
///
///   root -b -q 'be2_summary_C14.C()'
void be2_summary_C14()
{
   gStyle->SetOptStat(0);
   // REGENERATED 2026-08-28 on the d5 distributions (5 deg bins, scale over 25-105 deg), with
   // each potential scanned against ITS OWN 7.012 curve -- the consistent test. Cache and
   // luminosity are the adopted ones (cat5_s013, L = 72.5); only the binning and window differ
   // from the previous 10 deg / 25-135 version. Perey is included because it is adopted, and CH89
   // because it gives the best shape.
   const int NP = 5;
   const char *nm[NP] = {"KD03", "CH89", "Becchetti-Greenlees", "Perey", "Menet"};
   // the full name does not fit on the x axis and was clipped to "hetti-Greenlees"
   const char *nmShort[NP] = {"KD03", "CH89", "Becchetti-G.", "Perey", "Menet"};
   // B(E2)up [e2fm4] over f903 = 0, 0.1, 0.2, 0.3, 0.5, from f903_scan_C14.C on d5f*
   double B[NP][5] = {{40.5, 37.8, 35.3, 33.3, 29.4},
                      {49.0, 45.7, 42.8, 40.3, 35.6},
                      {43.1, 40.4, 37.6, 35.1, 31.3},
                      {44.2, 41.3, 38.5, 35.9, 32.1},
                      {43.7, 40.9, 38.1, 35.8, 31.8}};
   double rms[NP] = {0.290, 0.268, 0.278, 0.266, 0.302};      // best 7.012 shape over the scan
   // dip position from elastic_dip_C14.C; eps = L(dip)/L_scaler with L_scaler = 96.8 (clock)
   double dip[NP] = {63.0, 62.0, 58.0, 61.0, 58.0};
   double eps[NP] = {0.590, 0.805, 0.855, 0.736, 0.800};
   const double WU = 0.05940 * std::pow(14.0, 4.0 / 3.0);
   const double LIT = 5 * 1.8 * WU, LITe = 5 * 0.3 * WU;
   const double thDip = 57.95;
   int col[NP] = {kGray + 2, kAzure + 2, kGreen + 3, kRed + 1, kOrange + 8};
   const int ADOPTED = 3;        // Perey -- highlighted, it is the adopted potential

   auto *c = new TCanvas("cbe", "", 1400, 560); c->Divide(2, 1);
   c->cd(1); gPad->SetGridy(); gPad->SetBottomMargin(0.20);
   auto *fr = gPad->DrawFrame(-0.5, 0, NP - 0.5, 62);
   fr->SetTitle("B(E2) of the 7.012 level, 0^{+} #rightarrow 2^{+};;B(E2)#uparrow [e^{2}fm^{4}]");
   for (int i = 0; i < NP; ++i) fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin((double)i), nmShort[i]);
   fr->GetXaxis()->SetLabelSize(0.042);
   auto *bx = new TBox(-0.5, LIT - LITe, NP - 0.5, LIT + LITe);
   bx->SetFillColorAlpha(kRed + 1, 0.25); bx->Draw();
   auto *ll = new TLine(-0.5, LIT, NP - 0.5, LIT); ll->SetLineColor(kRed + 1); ll->SetLineWidth(3); ll->Draw();
   for (int i = 0; i < NP; ++i) {
      double hi = B[i][0], lo = B[i][4], cen = 0.5 * (hi + lo), half = 0.5 * (hi - lo);
      double tot = std::sqrt(half * half + std::pow(0.13 * cen, 2));
      auto *g = new TGraphErrors(1); g->SetPoint(0, i, cen); g->SetPointError(0, 0, tot);
      g->SetMarkerStyle(20); g->SetMarkerSize(2.2); g->SetMarkerColor(col[i]);
      g->SetLineColor(col[i]); g->SetLineWidth(3); g->Draw("P same");
      // the f903 band itself, as a bar
      auto *b2 = new TBox(i - 0.12, lo, i + 0.12, hi);
      b2->SetFillColorAlpha(col[i], 0.25); b2->SetLineColor(col[i]); b2->Draw("l");
      if (i == ADOPTED) { g->SetMarkerStyle(29); g->SetMarkerSize(3.4); }   // adopted: star
      TLatex t; t.SetTextSize(0.033); t.SetTextColor(col[i]); t.SetTextAlign(21);
      t.DrawLatex(i, cen + tot + 1.8, Form("%.0f #pm %.0f", cen, tot));
   }
   TLatex tr; tr.SetNDC(); tr.SetTextSize(0.040); tr.SetTextColor(kRed + 1);
   tr.DrawLatex(0.15, 0.27, "ENSDF 1.8(3) W.u. (decay) #times 5");
   tr.SetTextSize(0.032); tr.DrawLatex(0.15, 0.22, "star = adopted (Perey)");
   tr.SetTextColor(kBlack); tr.SetTextSize(0.034);
   tr.DrawLatex(0.15, 0.87, "bar = 6.903 blend band (f_{903} = 0 #rightarrow 0.5)");
   tr.DrawLatex(0.15, 0.82, "error bar = blend #oplus potential");

   // right: why BG, on the tests that do not involve B(E2) at all
   c->cd(2); gPad->SetGridx(); gPad->SetGridy();
   auto *f2 = gPad->DrawFrame(-1.5, 0.45, 7.5, 1.06);
   f2->SetTitle("the two tests with no DWBA normalisation in them;#theta_{dip}^{calc} - #theta_{dip}^{data} [deg];#varepsilon = L_{elastic} / L_{scaler}");
   for (int i = 0; i < NP; ++i) {
      auto *g = new TGraph(1); g->SetPoint(0, dip[i] - thDip, eps[i]);
      g->SetMarkerStyle(i == ADOPTED ? 29 : 20); g->SetMarkerSize(i == ADOPTED ? 3.6 : 2.4);
      g->SetMarkerColor(col[i]); g->Draw("P same");
      TLatex t; t.SetTextSize(0.031); t.SetTextColor(col[i]);
      // BG and Menet sit at the SAME dip offset (both 58.0 deg), so any parity-based staggering
      // puts their labels on top of each other. Place each one explicitly.
      const double LX[NP] = {-0.25,  0.30,  0.30,  0.30, 0.30};
      const double LY[NP] = { 0.030, 0.000, 0.014, 0.000, -0.016};
      const short  LA[NP] = { 32,    12,    12,    12,    12};
      t.SetTextAlign(LA[i]);
      t.DrawLatex(dip[i] - thDip + LX[i], eps[i] + LY[i], Form("%s (rms %.3f)", nm[i], rms[i]));
   }
   auto *v = new TLine(0, 0.45, 0, 1.06); v->SetLineColor(kRed + 1); v->SetLineStyle(2); v->SetLineWidth(2); v->Draw();
   auto *h1 = new TLine(-1.5, 1.0, 7.5, 1.0); h1->SetLineColor(kRed + 1); h1->SetLineWidth(2); h1->Draw();
   TLatex t2; t2.SetTextSize(0.034);
   t2.SetTextSize(0.031); t2.DrawLatex(-1.45, 1.015, "#varepsilon = 1: no unexplained loss");
   t2.DrawLatex(-1.45, 0.53, "dip correct and least unexplained");
   t2.DrawLatex(-1.45, 0.49, "inefficiency #rightarrow BG;");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "11_be2_summary.png");
   printf("\n  wrote %s11_be2_summary.png\n\n", out.Data());
}
