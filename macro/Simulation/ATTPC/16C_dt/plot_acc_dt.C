/// @file plot_acc_dt.C
/// @brief Draw the (d,t) acceptance from acceptance_dt.root. Separate from the macro that
/// computes it, so a plotting failure can never cost a 20-minute analysis again.
void plot_acc_dt(TString f = "data/acceptance_dt.root", TString png = "data/acceptance_dt.png")
{
   gStyle->SetOptStat(0);
   TFile *F = TFile::Open(f);
   if (!F || F->IsZombie()) { printf("cannot open %s\n", f.Data()); return; }
   const char *tg[5] = {"gs_s3001", "ex1_s3001", "ex2_s3001", "ex3_s3001", "ex4_s3001"};
   const double ex[5] = {0.0, 0.740, 3.103, 4.657, 6.358};
   const int col[5] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};
   auto *c = new TCanvas("c", "acceptance", 1100, 750);
   auto *leg = new TLegend(0.15, 0.15, 0.45, 0.37);
   leg->SetBorderSize(0);
   bool first = true;
   for (int i = 0; i < 5; ++i) {
      auto *h = (TH1D *)F->Get(TString("hAcc_") + tg[i]);
      if (!h) continue;
      h->SetDirectory(nullptr);
      h->SetLineColor(col[i]); h->SetMarkerColor(col[i]);
      h->SetMarkerStyle(20); h->SetMarkerSize(0.9); h->SetLineWidth(2);
      h->SetMinimum(0); h->SetMaximum(1.09);
      h->SetTitle("16C(d,t)15C acceptance, #theta_{lab} 8-56 deg;#theta_{cm} [deg];acceptance");
      h->Draw(first ? "E1" : "E1 SAME");
      leg->AddEntry(h, TString::Format("E_{x} = %.3f MeV", ex[i]), "lp");
      first = false;
   }
   leg->Draw();
   c->SaveAs(png);
   printf("plot -> %s\n", png.Data());
}
