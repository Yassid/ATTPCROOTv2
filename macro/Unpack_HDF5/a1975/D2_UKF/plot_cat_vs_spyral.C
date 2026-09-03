/// @file plot_cat_vs_spyral.C
/// @brief 15C excitation energy: CATIMA vs Spyral InterpSolver vs the current production.
/// Spyral has ~1/4 the tracks, so counts (left) and unit-area shapes (right) are both shown --
/// only the shape panel compares Spyral fairly.
void plot_cat_vs_spyral(TString out = "plots/cat_vs_spyral.png")
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   auto *fC = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fS = TFile::Open("/mnt/f/a1975/caches/dt_kin_spyral.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   auto *tC = (TNtuple *)fC->Get("pk");
   auto *tS = (TNtuple *)fS->Get("pk");
   auto *tO = (TNtuple *)fO->Get("pk");
   const char *cg = "chi2ndf<1e9";

   auto mk = [](const char *n) {
      auto *h = new TH1D(n, ";E_{x}(^{15}C)  (MeV);", 110, -3, 8);
      h->SetLineWidth(2); h->Sumw2(); return h; };
   auto *hO = mk("hO"); auto *hC = mk("hC"); auto *hS = mk("hS");
   tO->Draw("ex>>hO", cg, "goff"); tC->Draw("ex>>hC", cg, "goff"); tS->Draw("ex>>hS", "", "goff");
   hO->SetLineColor(kAzure + 2); hC->SetLineColor(kOrange + 8);
   hS->SetLineColor(kTeal + 3); hS->SetLineStyle(2);

   auto *c = new TCanvas("c", "", 1800, 700); c->Divide(2, 1);
   const double lv[] = {0.0, 0.740, 3.103, 4.220, 4.657};

   c->cd(1); gPad->SetLeftMargin(0.12);
   hO->SetTitle("counts (Spyral has ~1/4 the tracks)");
   hO->GetYaxis()->SetTitle("counts / 100 keV");
   hO->SetMaximum(1.25 * TMath::Max(hO->GetMaximum(), hC->GetMaximum())); hO->SetMinimum(0);
   hO->Draw("E0"); hC->Draw("E0 same"); hS->Draw("E0 same");
   for (double x : lv) { auto *l = new TLine(x,0,x,hO->GetMaximum());
      l->SetLineStyle(3); l->SetLineColor(kGray+2); l->Draw(); }
   auto *l1 = new TLegend(0.55,0.70,0.89,0.88); l1->SetBorderSize(0); l1->SetFillStyle(0);
   l1->SetTextSize(0.030);
   l1->AddEntry(hO, Form("matFX off  (%.0f)", hO->GetEntries()), "l");
   l1->AddEntry(hC, Form("CATIMA     (%.0f)", hC->GetEntries()), "l");
   l1->AddEntry(hS, Form("Spyral     (%.0f)", hS->GetEntries()), "l");
   l1->Draw();

   c->cd(2); gPad->SetLeftMargin(0.12);
   auto *nO=(TH1D*)hO->Clone("nO"); auto *nC=(TH1D*)hC->Clone("nC"); auto *nS=(TH1D*)hS->Clone("nS");
   for (auto *h : {nO,nC,nS}) if (h->Integral()>0) h->Scale(1.0/h->Integral());
   nO->SetTitle("unit area - the fair comparison");
   nO->GetYaxis()->SetTitle("fraction / 100 keV");
   nO->SetMaximum(1.25 * TMath::Max(nS->GetMaximum(), TMath::Max(nO->GetMaximum(), nC->GetMaximum())));
   nO->SetMinimum(0);
   nO->Draw("hist"); nC->Draw("hist same"); nS->Draw("hist same");
   for (double x : lv) { auto *l = new TLine(x,0,x,nO->GetMaximum());
      l->SetLineStyle(3); l->SetLineColor(kGray+2); l->Draw();
      auto *t=new TLatex(x+0.04, 1.14*nO->GetMaximum()/1.25, Form("%.3f",x));
      t->SetTextSize(0.024); t->SetTextColor(kGray+3); t->SetTextAngle(90); t->Draw(); }
   c->SaveAs(out);
   printf("\nwrote %s\n", out.Data());
   printf("%-8s %10s %10s %10s\n","level","matFXoff","CATIMA","Spyral");
   for (double x : lv) {
      int b1=hO->FindBin(x-0.25), b2=hO->FindBin(x+0.25);
      printf("%-8.3f %10.1f %10.1f %10.1f   (per 1000 tracks)\n", x,
             1000*hO->Integral(b1,b2)/hO->GetEntries(),
             1000*hC->Integral(b1,b2)/hC->GetEntries(),
             1000*hS->Integral(b1,b2)/hS->GetEntries());
   }
}
