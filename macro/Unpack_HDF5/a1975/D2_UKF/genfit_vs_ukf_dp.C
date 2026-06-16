/// @file genfit_vs_ukf_dp.C
/// @brief Overlay the 16C(d,p)17C proton spectrum: GENFIT vs UKF, same runs, same
/// proton gate. Reads two cached ntuples (pk) produced by ex_dp_a1975.C — one from
/// the genfit fits, one from the UKF fits (fitSuffix="_UKF"). Each is peak-normalized;
/// Ex, theta_lab and KE-vs-theta are overlaid. Saves plots/genfit_vs_ukf_dp.png.
///
///   root -b -q 'genfit_vs_ukf_dp.C("proton_kin_gf10.root","proton_kin_ukf10.root")'

void genfit_vs_ukf_dp(const char *gfFile = "proton_kin_gf10.root", const char *ukfFile = "proton_kin_ukf10.root")
{
   gStyle->SetOptStat(0);
   TFile *fg = TFile::Open(gfFile);
   TFile *fu = TFile::Open(ukfFile);
   TNtuple *tg = fg ? (TNtuple *)fg->Get("pk") : nullptr;
   TNtuple *tu = fu ? (TNtuple *)fu->Get("pk") : nullptr;
   if (!tg || !tu) {
      printf("missing ntuple (gf=%p ukf=%p)\n", (void *)tg, (void *)tu);
      return;
   }

   TH1F *exg = new TH1F("exg", "^{17}C E_{x}: genfit vs UKF;E_{x}(^{17}C) [MeV];norm.", 100, -10, 20);
   TH1F *exu = new TH1F("exu", "", 100, -10, 20);
   TH1F *thg = new TH1F("thg", "proton #theta_{lab};#theta_{lab} [deg];norm.", 60, 0, 180);
   TH1F *thu = new TH1F("thu", "", 60, 0, 180);
   tg->Draw("ex>>exg", "", "goff");
   tu->Draw("ex>>exu", "", "goff");
   tg->Draw("theta>>thg", "", "goff");
   tu->Draw("theta>>thu", "", "goff");
   for (auto h : {exg, thg}) {
      h->SetLineColor(kBlue + 1);
      h->SetLineWidth(2);
   }
   for (auto h : {exu, thu}) {
      h->SetLineColor(kRed);
      h->SetLineWidth(2);
   }
   for (auto h : {exg, exu, thg, thu})
      if (h->Integral() > 0)
         h->Scale(1. / h->Integral());

   long ng = tg->GetEntries(), nu = tu->GetEntries();
   double bg = tg->GetEntries("theta>90") * 100.0 / std::max(1L, ng);
   double bu = tu->GetEntries("theta>90") * 100.0 / std::max(1L, nu);
   printf("genfit candidates=%ld (backward %.1f%%)   UKF candidates=%ld (backward %.1f%%)\n", ng, bg, nu, bu);

   TCanvas *c = new TCanvas("c", "gf_vs_ukf", 1300, 560);
   c->Divide(2, 1);
   c->cd(1);
   exg->SetMaximum(1.25 * std::max(exg->GetMaximum(), exu->GetMaximum()));
   exg->Draw("hist");
   exu->Draw("hist same");
   TLegend *l = new TLegend(0.6, 0.7, 0.88, 0.88);
   l->AddEntry(exg, Form("genfit (%ld)", ng), "l");
   l->AddEntry(exu, Form("UKF (%ld)", nu), "l");
   l->Draw();
   TLine *l0 = new TLine(0, 0, 0, exg->GetMaximum());
   l0->SetLineStyle(2);
   l0->SetLineColor(kGray + 2);
   l0->Draw();
   c->cd(2);
   thg->SetMaximum(1.25 * std::max(thg->GetMaximum(), thu->GetMaximum()));
   thg->Draw("hist");
   thu->Draw("hist same");
   TLine *l90 = new TLine(90, 0, 90, thg->GetMaximum());
   l90->SetLineStyle(2);
   l90->SetLineColor(kGray + 2);
   l90->Draw();
   c->SaveAs("plots/genfit_vs_ukf_dp.png");
   printf("saved plots/genfit_vs_ukf_dp.png\n");
}
