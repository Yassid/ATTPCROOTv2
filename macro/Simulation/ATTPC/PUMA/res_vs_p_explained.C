/// @file res_vs_p_explained.C
/// @brief Momentum resolution vs |p| (UKF, GENFIT) with the two physics regimes
///        annotated: multiple scattering at low p, near-straight tracks at high p.
///        The U-shape has its minimum in the PUMA physics range. Values from the
///        momentum scan (dlc_ukf_genfit.C).
/// Run: root -b -q 'res_vs_p_explained.C'
void res_vs_p_explained(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/res_vs_p.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const int N = 3;
   double p[N]  = {150, 375, 600};
   double sG[N] = {23.0, 6.9, 9.3};   // GENFIT sigma [%]
   double sU[N] = {24.1, 15.3, 31.6}; // UKF sigma [%]

   auto *c = new TCanvas("c", "", 880, 600);
   c->SetLeftMargin(0.12); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetGrid();
   auto *gG = new TGraph(N, p, sG); auto *gU = new TGraph(N, p, sU);
   gG->SetTitle(";truth |p| [MeV/c];momentum resolution #sigma_{IQR}(p)/p  [%]");
   gG->SetMarkerStyle(21); gG->SetMarkerSize(1.7); gG->SetMarkerColor(kRed + 1); gG->SetLineColor(kRed + 1); gG->SetLineWidth(3); gG->SetLineStyle(2);
   gG->GetYaxis()->SetRangeUser(0, 40); gG->GetXaxis()->SetLimits(80, 680); gG->Draw("ALP");
   gU->SetMarkerStyle(20); gU->SetMarkerSize(1.7); gU->SetMarkerColor(kAzure + 2); gU->SetLineColor(kAzure + 2); gU->SetLineWidth(3); gU->Draw("LP");

   // regime shading
   auto *bxL = new TBox(80, 0, 250, 40); bxL->SetFillColorAlpha(kOrange - 9, 0.4); bxL->SetLineWidth(0); bxL->Draw();
   auto *bxH = new TBox(500, 0, 680, 40); bxH->SetFillColorAlpha(kAzure - 9, 0.4); bxH->SetLineWidth(0); bxH->Draw();
   gG->Draw("LP"); gU->Draw("LP"); gPad->RedrawAxis();

   auto *lg = new TLegend(0.40, 0.74, 0.72, 0.88); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(gU, "UKF", "lp"); lg->AddEntry(gG, "GENFIT", "lp"); lg->Draw();
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.030);
   tx->SetTextColor(kOrange + 8); tx->DrawLatex(112, 34, "multiple");
   tx->DrawLatex(112, 31.5, "scattering");
   tx->SetTextColor(kAzure + 3); tx->DrawLatex(520, 34, "near-straight");
   tx->DrawLatex(520, 31.5, "tracks");
   tx->SetTextColor(kGreen + 3); tx->SetTextSize(0.032); tx->DrawLatex(300, 4.0, "measurement-limited");
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.043);
   tt->DrawLatex(0.12, 0.93, "Momentum resolution vs |p| --- the U-shape");
   c->SaveAs(out);
   printf("RES_VS_P_DONE\n");
}
