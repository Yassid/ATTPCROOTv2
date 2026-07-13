/// @file pra_algo_plot.C
/// @brief Bar chart comparing PUMA track-finding algorithms at 375 MeV/c on the
///        2-pion channel: HC spatial clustering vs Riemann circle fit vs standalone
///        RANSAC 2D-circle. Global circle-consensus methods merge the back-to-back
///        pair; only clustering separates them. Values from pra_efficiency.C.
/// Run: root -b -q 'pra_algo_plot.C'
void pra_algo_plot(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pra_algo_cmp.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickY(1);
   const int NA = 3;
   const char *alg[NA] = {"HC clustering", "Riemann", "RANSAC circle"};
   double eff[NA]  = {99.3, 48.0, 51.5};
   double good[NA] = {86, 8, 11};
   double pur[NA]  = {99.2, 61.9, 55.4};

   auto *hEff = new TH1F("hEff", ";;fraction [%]", NA, 0, NA);
   auto *hGood = new TH1F("hGood", "", NA, 0, NA);
   auto *hPur = new TH1F("hPur", "", NA, 0, NA);
   for (int i = 0; i < NA; ++i) { hEff->SetBinContent(i + 1, eff[i]); hGood->SetBinContent(i + 1, good[i]); hPur->SetBinContent(i + 1, pur[i]);
      hEff->GetXaxis()->SetBinLabel(i + 1, alg[i]); }
   hEff->GetXaxis()->SetLabelSize(0.05);
   hEff->SetBarWidth(0.26); hEff->SetBarOffset(0.08); hEff->SetFillColor(kAzure + 2);
   hGood->SetBarWidth(0.26); hGood->SetBarOffset(0.37); hGood->SetFillColor(kGreen + 2);
   hPur->SetBarWidth(0.26); hPur->SetBarOffset(0.66); hPur->SetFillColor(kOrange + 7);

   auto *c = new TCanvas("c", "", 850, 600); c->SetLeftMargin(0.11); c->SetBottomMargin(0.10); c->SetTopMargin(0.09);
   hEff->SetMaximum(115); hEff->SetMinimum(0);
   hEff->Draw("bar"); hGood->Draw("bar same"); hPur->Draw("bar same");
   auto *lg = new TLegend(0.5, 0.75, 0.88, 0.90); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hEff, "track-finding efficiency", "f");
   lg->AddEntry(hPur, "purity", "f");
   lg->AddEntry(hGood, "clean 2-track events", "f");
   lg->Draw();
   // value labels
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.028); tx->SetTextAlign(21);
   for (int i = 0; i < NA; ++i) {
      tx->SetTextColor(kAzure + 2); tx->DrawLatex(i + 0.21, eff[i] + 2, Form("%.0f", eff[i]));
      tx->SetTextColor(kOrange + 8); tx->DrawLatex(i + 0.50, pur[i] + 2, Form("%.0f", pur[i]));
      tx->SetTextColor(kGreen + 3); tx->DrawLatex(i + 0.79, good[i] + 2, Form("%.0f", good[i]));
   }
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.042);
   tt->DrawLatex(0.11, 0.93, "PUMA track-finding algorithms @ 375 MeV/c (2 pions)");
   auto *tn = new TLatex(); tn->SetNDC(); tn->SetTextFont(62); tn->SetTextSize(0.03); tn->SetTextColor(kRed + 1);
   tn->DrawLatex(0.40, 0.20, "global circle fits merge the pair");
   c->SaveAs(out);
   printf("PRAALGO_PLOT_DONE\n");
}
