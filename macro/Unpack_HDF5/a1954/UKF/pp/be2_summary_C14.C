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
   const int NP = 3;
   const char *nm[NP] = {"KD03", "Menet", "Becchetti-Greenlees"};
   // B(E2)up over f903 = 0, 0.1, 0.2, 0.3, 0.5, from f903_scan_C14.C
   double B[NP][5] = {{45.7, 43.9, 40.9, 38.2, 32.7},
                      {37.1, 35.8, 33.4, 31.2, 26.6},
                      {34.6, 33.2, 30.9, 28.9, 24.8}};
   double rms[NP] = {0.449, 0.395, 0.365};      // best 7.012 shape over the scan
   double dip[NP] = {63.0, 58.0, 58.0}, eps[NP] = {0.59, 0.80, 0.86};
   const double WU = 0.05940 * std::pow(14.0, 4.0 / 3.0);
   const double LIT = 5 * 1.8 * WU, LITe = 5 * 0.3 * WU;
   const double thDip = 57.95;
   int col[NP] = {kGray + 2, kAzure + 2, kGreen + 3};

   auto *c = new TCanvas("cbe", "", 1400, 560); c->Divide(2, 1);
   c->cd(1); gPad->SetGridy(); gPad->SetBottomMargin(0.20);
   auto *fr = gPad->DrawFrame(-0.5, 0, NP - 0.5, 60);
   fr->SetTitle("B(E2; 0^{+}#rightarrow2^{+}) of the 7.012 level;;B(E2)#uparrow [e^{2}fm^{4}]");
   for (int i = 0; i < NP; ++i) fr->GetXaxis()->SetBinLabel(fr->GetXaxis()->FindBin((double)i), nm[i]);
   fr->GetXaxis()->SetLabelSize(0.050);
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
      TLatex t; t.SetTextSize(0.036); t.SetTextColor(col[i]); t.SetTextAlign(21);
      t.DrawLatex(i, cen + tot + 2.5, Form("%.0f #pm %.0f", cen, tot));
   }
   TLatex tr; tr.SetNDC(); tr.SetTextSize(0.040); tr.SetTextColor(kRed + 1);
   tr.DrawLatex(0.15, 0.22, "ENSDF 1.8(3) W.u. (decay) #times 5");
   tr.SetTextColor(kBlack); tr.SetTextSize(0.034);
   tr.DrawLatex(0.15, 0.87, "bar = 6.903 blend band (f_{903} = 0 #rightarrow 0.5)");
   tr.DrawLatex(0.15, 0.82, "error bar = blend #oplus potential");

   // right: why BG, on the tests that do not involve B(E2) at all
   c->cd(2); gPad->SetGridx(); gPad->SetGridy();
   auto *f2 = gPad->DrawFrame(-1.0, 0.5, 6.5, 1.05);
   f2->SetTitle("the two tests with no DWBA normalisation in them;#theta_{dip}^{calc} - #theta_{dip}^{data} [deg];#varepsilon = L_{elastic} / L_{scaler}");
   for (int i = 0; i < NP; ++i) {
      auto *g = new TGraph(1); g->SetPoint(0, dip[i] - thDip, eps[i]);
      g->SetMarkerStyle(20); g->SetMarkerSize(2.4); g->SetMarkerColor(col[i]); g->Draw("P same");
      TLatex t; t.SetTextSize(0.038); t.SetTextColor(col[i]);
      t.DrawLatex(dip[i] - thDip + 0.2, eps[i] - 0.005, Form("%s  (7.012 rms %.3f)", nm[i], rms[i]));
   }
   auto *v = new TLine(0, 0.5, 0, 1.05); v->SetLineColor(kRed + 1); v->SetLineStyle(2); v->SetLineWidth(2); v->Draw();
   auto *h1 = new TLine(-1.0, 1.0, 6.5, 1.0); h1->SetLineColor(kRed + 1); h1->SetLineWidth(2); h1->Draw();
   TLatex t2; t2.SetTextSize(0.034);
   t2.DrawLatex(-0.85, 1.015, "#varepsilon = 1: no unexplained loss");
   t2.DrawLatex(1.2, 0.56, "dip correct #rightarrow");
   t2.DrawLatex(1.2, 0.53, "and least unexplained inefficiency #rightarrow BG");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "11_be2_summary.png");
   printf("\n  wrote %s11_be2_summary.png\n\n", out.Data());
}
