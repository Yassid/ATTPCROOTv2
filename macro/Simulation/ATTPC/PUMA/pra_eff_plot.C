/// @file pra_eff_plot.C
/// @brief Plot PUMA PRA (smooth-3D clustering) track-finding efficiency and
///        clean-2-track ("good") event fraction vs momentum, showing the two
///        failure modes: over-segmentation (split) at low p, track merging at
///        high p. Values from pra_efficiency.C on the momentum scan.
/// Run: root -b -q 'pra_eff_plot.C'
void pra_eff_plot(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pra_eff_vs_p.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const int N = 6;
   double p[N]   = {150, 250, 375, 500, 700, 900};
   double eff[N] = {99.3, 97.3, 99.0, 97.4, 74.0, 62.2};
   double good[N]= {5, 9, 96, 95, 48, 24};
   double split[N]={94, 86, 2, 0, 0, 0};
   double merge[N]={0, 0, 2, 5, 52, 76};

   auto *c = new TCanvas("c", "", 900, 620);
   c->SetLeftMargin(0.12); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetGrid();
   auto *gEff = new TGraph(N, p, eff);
   auto *gGood = new TGraph(N, p, good);
   auto *gSplit = new TGraph(N, p, split);
   auto *gMerge = new TGraph(N, p, merge);
   gEff->SetTitle(";truth |p| [MeV/c];fraction [%]");
   gEff->SetMarkerStyle(20); gEff->SetMarkerSize(1.4); gEff->SetMarkerColor(kAzure + 2); gEff->SetLineColor(kAzure + 2); gEff->SetLineWidth(3);
   gEff->GetYaxis()->SetRangeUser(0, 108); gEff->GetXaxis()->SetLimits(100, 950);
   gEff->Draw("ALP");
   gGood->SetMarkerStyle(21); gGood->SetMarkerSize(1.4); gGood->SetMarkerColor(kGreen + 2); gGood->SetLineColor(kGreen + 2); gGood->SetLineWidth(3); gGood->Draw("LP");
   gSplit->SetMarkerStyle(24); gSplit->SetMarkerSize(1.2); gSplit->SetMarkerColor(kOrange + 7); gSplit->SetLineColor(kOrange + 7); gSplit->SetLineWidth(2); gSplit->SetLineStyle(2); gSplit->Draw("LP");
   gMerge->SetMarkerStyle(25); gMerge->SetMarkerSize(1.2); gMerge->SetMarkerColor(kRed + 1); gMerge->SetLineColor(kRed + 1); gMerge->SetLineWidth(2); gMerge->SetLineStyle(2); gMerge->Draw("LP");

   auto *lg = new TLegend(0.44, 0.60, 0.80, 0.88); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(gEff, "track-finding efficiency", "lp");
   lg->AddEntry(gGood, "clean 2-track (\"good\") events", "lp");
   lg->AddEntry(gSplit, "split (>2 tracks)", "lp");
   lg->AddEntry(gMerge, "merged (1 track)", "lp");
   lg->Draw();
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.032);
   tx->SetTextColor(kOrange + 8); tx->DrawLatex(150, 60, "#leftarrow over-");
   tx->DrawLatex(150, 54, "segmentation");
   tx->SetTextColor(kRed + 1); tx->DrawLatex(760, 60, "merging #rightarrow");
   tx->SetTextColor(kBlack); tx->SetTextSize(0.036); tx->DrawLatex(360, 101, "sweet spot");
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.045);
   tt->DrawLatex(0.13, 0.93, "PUMA PRA (smooth-3D) efficiency vs momentum");
   c->SaveAs(out);
   printf("PRAEFF_PLOT_DONE\n");
}
