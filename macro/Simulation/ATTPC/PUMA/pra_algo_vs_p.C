/// @file pra_algo_vs_p.C
/// @brief Clean-2-track ("good") event fraction vs momentum for the two clustering
///        track finders on the PUMA 2-pion channel (same reco config): TriplClust
///        (AtTrackFinderTC, default t) vs HDBSCAN (AtTrackFinderHDBSCAN). HDBSCAN
///        stays robust at high p where TriplClust-default merges; a star shows
///        TriplClust recovers at 900 with a lowered cluster distance (t=2).
/// Run: root -b -q 'pra_algo_vs_p.C'
void pra_algo_vs_p(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pra_algo_vs_p.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const int N = 6;
   double p[N]   = {150, 250, 375, 500, 700, 900};
   double tc[N]  = {1, 6, 86, 88, 50, 21};   // TriplClust (default t)
   double hd[N]  = {26, 12, 86, 90, 89, 83}; // HDBSCAN

   auto *c = new TCanvas("c", "", 900, 620);
   c->SetLeftMargin(0.11); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetGrid();
   auto *gTC = new TGraph(N, p, tc);
   auto *gHD = new TGraph(N, p, hd);
   gTC->SetTitle(";truth |p| [MeV/c];clean 2-track events [%]");
   gTC->SetMarkerStyle(21); gTC->SetMarkerSize(1.5); gTC->SetMarkerColor(kGreen + 2); gTC->SetLineColor(kGreen + 2); gTC->SetLineWidth(3); gTC->SetLineStyle(2);
   gTC->GetYaxis()->SetRangeUser(0, 108); gTC->GetXaxis()->SetLimits(100, 950);
   gTC->Draw("ALP");
   gHD->SetMarkerStyle(20); gHD->SetMarkerSize(1.5); gHD->SetMarkerColor(kAzure + 2); gHD->SetLineColor(kAzure + 2); gHD->SetLineWidth(3);
   gHD->Draw("LP");
   // TriplClust tuned (t=2) recovery point at 900
   auto *gTune = new TGraph(1); gTune->SetPoint(0, 900, 82);
   gTune->SetMarkerStyle(29); gTune->SetMarkerSize(2.6); gTune->SetMarkerColor(kOrange + 8); gTune->Draw("P");

   // kaon region band (~600 MeV/c)
   auto *bx = new TBox(560, 0, 640, 108); bx->SetFillColorAlpha(kGray, 0.25); bx->SetLineWidth(0); bx->Draw();
   gTC->Draw("LP"); gHD->Draw("LP"); gTune->Draw("P"); gPad->RedrawAxis();

   auto *lg = new TLegend(0.36, 0.18, 0.74, 0.40); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(gHD, "HDBSCAN (density)", "lp");
   lg->AddEntry(gTC, "TriplClust (default t)", "lp");
   lg->AddEntry(gTune, "TriplClust tuned (t=2) @900", "p");
   lg->Draw();
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.030); tx->SetTextColor(kGray + 3); tx->SetTextAngle(90);
   tx->DrawLatex(600, 44, "PUMA kaons");
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.044);
   tt->DrawLatex(0.11, 0.93, "PUMA track-finding: HDBSCAN vs TriplClust vs |p|");
   auto *tn = new TLatex(); tn->SetNDC(); tn->SetTextFont(62); tn->SetTextSize(0.03); tn->SetTextColor(kRed + 1);
   tn->DrawLatex(0.58, 0.55, "TriplClust-default");
   tn->DrawLatex(0.58, 0.51, "merges near-straight");
   c->SaveAs(out);
   printf("ALGO_VS_P_DONE\n");
}
