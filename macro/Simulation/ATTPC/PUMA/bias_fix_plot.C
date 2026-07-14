/// @file bias_fix_plot.C
/// @brief Momentum bias vs |p| before/after the deterministic vertex material
///        correction (Cu trap + Al cryostat energy loss). Values from the
///        momentum scan (compare_ukf_genfit_test8.C / scan_momentum_plot.C).
/// Run: root -b -q 'bias_fix_plot.C'
void bias_fix_plot(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/bias_fix.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const int N = 6;
   double p[N]    = {150, 250, 375, 500, 700, 900};
   double raw[N]  = {-20.3, -10.9, -4.0, -8.2, -2.6, -5.0}; // GENFIT median bias, no correction (900 UKF diverges; GENFIT shown)
   double corr[N] = {0.3, -2.2, -1.5, -0.5, 0.9, -2.0};      // with material correction (GENFIT)

   auto *c = new TCanvas("c", "", 880, 600);
   c->SetLeftMargin(0.12); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetGrid();
   auto *gR = new TGraph(N, p, raw);
   auto *gC = new TGraph(N, p, corr);
   gR->SetTitle(";truth |p| [MeV/c];momentum bias  (p_{fit}-p_{truth})/p_{truth} [%]");
   gR->SetMarkerStyle(24); gR->SetMarkerSize(1.5); gR->SetMarkerColor(kRed + 1); gR->SetLineColor(kRed + 1); gR->SetLineWidth(3); gR->SetLineStyle(2);
   gR->GetYaxis()->SetRangeUser(-24, 8); gR->GetXaxis()->SetLimits(100, 950);
   gR->Draw("ALP");
   gC->SetMarkerStyle(20); gC->SetMarkerSize(1.5); gC->SetMarkerColor(kAzure + 2); gC->SetLineColor(kAzure + 2); gC->SetLineWidth(3);
   gC->Draw("LP");
   auto *l0 = new TLine(100, 0, 950, 0); l0->SetLineColor(kGray + 2); l0->SetLineStyle(3); l0->Draw();
   // +/-2% band
   auto *bx = new TBox(100, -2, 950, 2); bx->SetFillColorAlpha(kGreen - 9, 0.35); bx->SetLineWidth(0); bx->Draw();
   gR->Draw("LP"); gC->Draw("LP"); gPad->RedrawAxis();

   auto *lg = new TLegend(0.40, 0.20, 0.78, 0.38); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(gR, "before (in-gas momentum)", "lp");
   lg->AddEntry(gC, "after material correction", "lp");
   lg->Draw();
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.043);
   tt->DrawLatex(0.12, 0.93, "Momentum bias fixed by the Cu-trap/cryostat correction");
   auto *tn = new TLatex(); tn->SetTextFont(62); tn->SetTextSize(0.032); tn->SetTextColor(kGreen + 3);
   tn->DrawLatex(620, 3.0, "#pm2% band");
   c->SaveAs(out);
   printf("BIASFIX_PLOT_DONE\n");
}
