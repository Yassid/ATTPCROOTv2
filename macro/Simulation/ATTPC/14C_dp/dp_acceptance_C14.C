/// @file dp_acceptance_C14.C
/// @brief Acceptance of 14C(d,p)15C vs theta_cm, one curve per configuration.
///
///   root -b -q 'dp_acceptance_C14.C()'
///
/// The overall acceptance of this channel is 0.63-0.71 in every cell of the matrix, which makes
/// the field look nearly free. It is not: the loss is concentrated at small theta_cm -- the
/// backward-going, low-energy proton -- which is exactly where a transfer angular distribution has
/// its yield. Plotting the ratio to the present configuration (right panel) is what makes the cost
/// visible; the left panel alone would let 7 T pass as harmless.

void dp_acceptance_C14(TString root = "/mnt/f/a1954_C14dp_hf", TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const int NC = 6;
   const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
   const char *CLB[NC] = {"2.85 T, AT-TPC", "2.85 T, 2 mm", "4 T, AT-TPC", "4 T, 2 mm", "7 T, AT-TPC", "7 T, 2 mm"};
   const int COL[NC] = {kGray + 3, kGray + 1, kAzure + 2, kAzure - 4, kRed + 1, kOrange + 7};
   const int MK[NC] = {20, 24, 21, 25, 22, 26};

   // rebin the 5-degree acceptance histograms into wider bins, so the backward end is not noise
   const int NB = 9;
   const double lo[NB] = {5, 10, 20, 30, 45, 60, 80, 110, 140};
   const double hi[NB] = {10, 20, 30, 45, 60, 80, 110, 140, 180};

   double acc[NC][NB], err[NC][NB];
   for (int c = 0; c < NC; ++c)
      for (int b = 0; b < NB; ++b) acc[c][b] = err[c][b] = -1;

   for (int c = 0; c < NC; ++c) {
      TString f = gSystem->GetFromPipe(
         TString::Format("ls %s/%s/acceptance_gs_s*_%s.root 2>/dev/null | head -1", root.Data(), CFG[c], CFG[c]));
      f = f.Strip(TString::kBoth);
      if (f.IsNull()) continue;
      TFile *fa = TFile::Open(f);
      if (!fa || fa->IsZombie()) continue;
      TH1D *g = nullptr, *r = nullptr;
      TIter nx(fa->GetListOfKeys());
      while (auto *k = (TKey *)nx()) {
         TString n = k->GetName();
         if (n.BeginsWith("hGen_")) g = (TH1D *)fa->Get(n);
         if (n.BeginsWith("hRec_")) r = (TH1D *)fa->Get(n);
      }
      if (!g || !r) { fa->Close(); continue; }
      for (int b = 0; b < NB; ++b) {
         double gg = g->Integral(g->FindBin(lo[b] + 0.01), g->FindBin(hi[b] - 0.01));
         double rr = r->Integral(r->FindBin(lo[b] + 0.01), r->FindBin(hi[b] - 0.01));
         if (gg <= 0) continue;
         acc[c][b] = rr / gg;
         err[c][b] = std::sqrt(std::max(1e-9, acc[c][b] * (1 - acc[c][b]) / gg)); // binomial
      }
      fa->Close();
   }

   auto *cv = new TCanvas("dpacc", "dpacc", 1500, 620);
   cv->Divide(2, 1);

   for (int pad = 1; pad <= 2; ++pad) {
      cv->cd(pad);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.14);
      gPad->SetGridy();
      auto *fr = new TH1F(Form("fa%d", pad), pad == 1 ? ";#theta_{cm} [deg];acceptance"
                                                      : ";#theta_{cm} [deg];acceptance / present configuration",
                          100, 0, 180);
      fr->GetYaxis()->SetRangeUser(0, pad == 1 ? 1.05 : 1.35);
      fr->GetYaxis()->SetTitleSize(0.045);
      fr->GetYaxis()->SetTitleOffset(1.25);
      fr->GetXaxis()->SetTitleSize(0.045);
      fr->SetLineColor(kWhite);
      fr->Draw();
      auto *lg = new TLegend(pad == 1 ? 0.55 : 0.16, 0.17, pad == 1 ? 0.90 : 0.52, 0.42);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      lg->SetTextSize(0.036);
      for (int c = 0; c < NC; ++c) {
         auto *gph = new TGraphErrors();
         for (int b = 0; b < NB; ++b) {
            if (acc[c][b] < 0) continue;
            double y = acc[c][b], e = err[c][b];
            if (pad == 2) {
               if (acc[0][b] <= 0) continue;
               e = (y / acc[0][b]) * std::sqrt(std::pow(err[c][b] / std::max(1e-9, y), 2) +
                                               std::pow(err[0][b] / acc[0][b], 2));
               y /= acc[0][b];
            }
            int i = gph->GetN();
            gph->SetPoint(i, 0.5 * (lo[b] + hi[b]), y);
            gph->SetPointError(i, 0.5 * (hi[b] - lo[b]), e);
         }
         gph->SetLineColor(COL[c]);
         gph->SetMarkerColor(COL[c]);
         gph->SetMarkerStyle(MK[c]);
         gph->SetMarkerSize(1.3);
         gph->SetLineWidth(2);
         gph->Draw("pl same");
         if (pad == 1) lg->AddEntry(gph, CLB[c], "pl");
      }
      if (pad == 1) lg->Draw();
      if (pad == 2) {
         auto *l1 = new TLine(0, 1, 180, 1);
         l1->SetLineStyle(2);
         l1->SetLineColor(kGray + 2);
         l1->Draw();
      }
      // the transfer-peak band
      auto *bx = new TBox(8, 0, 30, pad == 1 ? 1.05 : 1.35);
      bx->SetFillColorAlpha(kGray + 1, 0.15);
      bx->SetLineWidth(0);
      bx->Draw();
      gPad->RedrawAxis();
      auto *tx = new TLatex(0.16, 0.90,
                            pad == 1 ? "#bf{A}  acceptance vs #theta_{cm}" : "#bf{B}  cost relative to 2.85 T, AT-TPC");
      tx->SetNDC();
      tx->SetTextSize(0.042);
      tx->Draw();
      auto *tb = new TLatex(0.235, pad == 1 ? 0.80 : 0.80, "#splitline{transfer}{peaks here}");
      tb->SetNDC();
      tb->SetTextSize(0.032);
      tb->SetTextColor(kGray + 3);
      tb->Draw();
   }
   cv->SaveAs(outDir + "dp_acceptance.png");
   printf("\n  wrote %sdp_acceptance.png\n\n", outDir.Data());
}
