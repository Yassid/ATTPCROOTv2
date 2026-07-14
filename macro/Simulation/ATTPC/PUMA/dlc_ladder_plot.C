/// @file dlc_ladder_plot.C
/// @brief DLC momentum recovery ladder @375 MeV/c (UKF): how each reconstruction
///        lever recovers the resolution that DLC charge dispersion destroys.
///        no-DLC reference vs DLC with max PSA -> +MultiFit PSA -> +HDBSCAN.
/// Run: root -b -q 'dlc_ladder_plot.C'
void dlc_ladder_plot(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/dlc_ladder.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickY(1);
   const int N = 5;
   const char *cfg[N] = {"no DLC (UKF)", "DLC max+ring (UKF)", "DLC MultiFit+ring (UKF)",
                         "DLC +HDBSCAN (UKF)", "DLC +HDBSCAN (GENFIT)"};
   double sig[N]  = {15.0, 18.9, 17.8, 15.1, 7.1};  // momentum sigma_IQR [%]
   double bias[N] = {0.1, 20.8, 8.9, 6.3, 0.1};     // |momentum bias| [%]

   auto *hS = new TH1F("hS", ";;value [%]", N, 0, N);
   auto *hB = new TH1F("hB", "", N, 0, N);
   for (int i = 0; i < N; ++i) { hS->SetBinContent(i + 1, sig[i]); hB->SetBinContent(i + 1, bias[i]); hS->GetXaxis()->SetBinLabel(i + 1, cfg[i]); }
   hS->GetXaxis()->SetLabelSize(0.028);
   hS->SetBarWidth(0.40); hS->SetBarOffset(0.06); hS->SetFillColor(kAzure + 2);
   hB->SetBarWidth(0.40); hB->SetBarOffset(0.52); hB->SetFillColor(kRed + 1);

   auto *c = new TCanvas("c", "", 900, 600); c->SetLeftMargin(0.10); c->SetBottomMargin(0.11); c->SetTopMargin(0.10);
   hS->SetMaximum(26); hS->SetMinimum(0);
   hS->Draw("bar"); hB->Draw("bar same");
   auto *lg = new TLegend(0.55, 0.78, 0.90, 0.89); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hS, "momentum resolution #sigma", "f"); lg->AddEntry(hB, "|momentum bias|", "f");
   lg->Draw();
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.028); tx->SetTextAlign(21);
   for (int i = 0; i < N; ++i) { tx->SetTextColor(kAzure + 2); tx->DrawLatex(i + 0.26, sig[i] + 0.6, Form("%.0f", sig[i]));
      tx->SetTextColor(kRed + 1); tx->DrawLatex(i + 0.72, bias[i] + 0.6, Form("%.0f", bias[i])); }
   auto *ref = new TLine(0, 15, N, 15); ref->SetLineColor(kAzure + 2); ref->SetLineStyle(2); ref->Draw();
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.042);
   tt->DrawLatex(0.10, 0.93, "Recovering (and beating) DLC momentum resolution @375 MeV/c");
   auto *tn = new TLatex(); tn->SetNDC(); tn->SetTextFont(62); tn->SetTextSize(0.027); tn->SetTextColor(kAzure + 3);
   tn->DrawLatex(0.30, 0.60, "no-DLC #sigma (15%)");
   tn->SetTextColor(kGreen + 3); tn->DrawLatex(0.63, 0.30, "GENFIT #sigma 7% = 2#times BETTER");
   tn->DrawLatex(0.63, 0.26, "than no-DLC (sub-pad gain)");
   c->SaveAs(out);
   printf("DLC_LADDER_DONE\n");
}
