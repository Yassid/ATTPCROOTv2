/// @file pra_kaon_cmp.C
/// @brief The decisive PUMA case: same-charge K+K+ at 600 MeV/c (near-straight,
///        mirror circles). TriplClust merges the pair; HDBSCAN separates it.
///        Values from pra_efficiency.C (targetPDG=321) on the 400-evt K+K+ sim.
/// Run: root -b -q 'pra_kaon_cmp.C'
void pra_kaon_cmp(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pra_kaon_cmp.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickY(1);
   const int NM = 4;
   const char *metric[NM] = {"efficiency", "purity", "clean 2-track", "merged"};
   double tc[NM] = {68.8, 77.2, 33, 57};
   double hd[NM] = {99.1, 98.8, 89, 0};

   auto *hTC = new TH1F("hTC", ";;fraction [%]", NM, 0, NM);
   auto *hHD = new TH1F("hHD", "", NM, 0, NM);
   for (int i = 0; i < NM; ++i) { hTC->SetBinContent(i + 1, tc[i]); hHD->SetBinContent(i + 1, hd[i]); hTC->GetXaxis()->SetBinLabel(i + 1, metric[i]); }
   hTC->GetXaxis()->SetLabelSize(0.055);
   hTC->SetBarWidth(0.40); hTC->SetBarOffset(0.06); hTC->SetFillColor(kGreen + 2);
   hHD->SetBarWidth(0.40); hHD->SetBarOffset(0.52); hHD->SetFillColor(kAzure + 2);

   auto *c = new TCanvas("c", "", 860, 600); c->SetLeftMargin(0.11); c->SetBottomMargin(0.10); c->SetTopMargin(0.11);
   hTC->SetMaximum(118); hTC->SetMinimum(0);
   hTC->Draw("bar"); hHD->Draw("bar same");
   auto *lg = new TLegend(0.50, 0.79, 0.88, 0.90); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hTC, "TriplClust", "f"); lg->AddEntry(hHD, "HDBSCAN", "f");
   lg->Draw();
   auto *tx = new TLatex(); tx->SetTextFont(62); tx->SetTextSize(0.030); tx->SetTextAlign(21);
   for (int i = 0; i < NM; ++i) {
      tx->SetTextColor(kGreen + 3); tx->DrawLatex(i + 0.26, tc[i] + 2, Form("%.0f", tc[i]));
      tx->SetTextColor(kAzure + 2); tx->DrawLatex(i + 0.72, hd[i] + 2, Form("%.0f", hd[i]));
   }
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.040);
   tt->DrawLatex(0.11, 0.94, "Same-charge K^{+}K^{+} @ 600 MeV/c track separation");
   auto *tn = new TLatex(); tn->SetNDC(); tn->SetTextFont(62); tn->SetTextSize(0.032); tn->SetTextColor(kRed + 1);
   tn->DrawLatex(0.58, 0.40, "TriplClust merges 57%");
   tn->DrawLatex(0.58, 0.36, "of the K^{+}K^{+} pairs");
   c->SaveAs(out);
   printf("KAON_CMP_DONE\n");
}
