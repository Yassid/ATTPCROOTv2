/// @file dlc_ukf_genfit.C
/// @brief UKF vs GENFIT momentum resolution & bias vs |p| with DLC ON, winning
///        config (HDBSCAN + MultiFit PSA + ring clustering). GENFIT exploits the
///        charge-weighted sub-pad DLC centroids (0.3 mm meas.) -> unbiased, sigma
///        7-9% in the PUMA physics range; UKF recovers the no-DLC sigma but keeps
///        a small bias and degrades on near-straight high-p tracks.
/// Run: root -b -q 'dlc_ukf_genfit.C'
void dlc_ukf_genfit(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/dlc_ukf_genfit.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const int N = 3;
   double p[N]    = {150, 375, 600};
   double sU[N]   = {22.4, 15.1, 42.1};  // UKF sigma
   double sG[N]   = {21.8, 7.1, 9.3};    // GENFIT sigma
   double bU[N]   = {-9.4, -6.3, -9.9};  // UKF bias
   double bG[N]   = {-2.8, -0.1, -0.4};  // GENFIT bias

   auto *c = new TCanvas("c", "", 1050, 500); c->Divide(2, 1);
   // sigma panel
   c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   auto *gU = new TGraph(N, p, sU), *gG = new TGraph(N, p, sG);
   gU->SetTitle("momentum resolution (DLC on);truth |p| [MeV/c];#sigma_{IQR}(p)/p  [%]");
   gU->SetMarkerStyle(20); gU->SetMarkerSize(1.6); gU->SetMarkerColor(kAzure + 2); gU->SetLineColor(kAzure + 2); gU->SetLineWidth(3);
   gU->GetYaxis()->SetRangeUser(0, 48); gU->GetXaxis()->SetLimits(100, 650); gU->Draw("ALP");
   gG->SetMarkerStyle(21); gG->SetMarkerSize(1.6); gG->SetMarkerColor(kRed + 1); gG->SetLineColor(kRed + 1); gG->SetLineWidth(3); gG->SetLineStyle(2); gG->Draw("LP");
   { auto *l = new TLine(100, 15, 650, 15); l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->Draw(); }
   { auto *lg = new TLegend(0.42, 0.72, 0.86, 0.88); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(gU, "UKF", "lp"); lg->AddEntry(gG, "GENFIT", "lp"); lg->Draw(); }
   { auto *t = new TLatex(); t->SetTextFont(62); t->SetTextSize(0.033); t->SetTextColor(kGray + 3); t->DrawLatex(400, 16.5, "no-DLC #sigma"); }
   // bias panel
   c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   auto *bUg = new TGraph(N, p, bU), *bGg = new TGraph(N, p, bG);
   bUg->SetTitle("momentum bias (DLC on);truth |p| [MeV/c];median bias [%]");
   bUg->SetMarkerStyle(20); bUg->SetMarkerSize(1.6); bUg->SetMarkerColor(kAzure + 2); bUg->SetLineColor(kAzure + 2); bUg->SetLineWidth(3);
   bUg->GetYaxis()->SetRangeUser(-14, 6); bUg->GetXaxis()->SetLimits(100, 650); bUg->Draw("ALP");
   bGg->SetMarkerStyle(21); bGg->SetMarkerSize(1.6); bGg->SetMarkerColor(kRed + 1); bGg->SetLineColor(kRed + 1); bGg->SetLineWidth(3); bGg->SetLineStyle(2); bGg->Draw("LP");
   { auto *l = new TLine(100, 0, 650, 0); l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->Draw(); }
   { auto *lg = new TLegend(0.42, 0.20, 0.86, 0.36); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(bUg, "UKF", "lp"); lg->AddEntry(bGg, "GENFIT (unbiased)", "lp"); lg->Draw(); }
   c->SaveAs(out);
   printf("DLC_UKFGF_DONE\n");
}
