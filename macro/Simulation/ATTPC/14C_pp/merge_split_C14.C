/// @file merge_split_C14.C
/// @brief Sum the per-seed acceptance_split_C14.C outputs and compare UKF against GENFIT.
///
/// Summing the RAW numerators and denominator before dividing (rather than averaging five
/// ratios) is what keeps the binomial errors honest at small bin populations.
///
///   root -b -q 'merge_split_C14.C("diagnostics/split","gs")'

void merge_split_C14(TString dir = "diagnostics/split", TString tag = "gs", Int_t s0 = 1001, Int_t s1 = 1005)
{
   gStyle->SetOptStat(0);
   const char *fitName[2] = {"ukf", "gf"};
   const char *fitLbl[2] = {"UKF", "GENFIT"};
   enum { CAT_NOTRK = 0, CAT_KEBAD, CAT_CHI2, CAT_BOTH, CAT_DTH, CAT_RATIO, CAT_GOOD, NCAT };
   const char *catName[NCAT] = {"NOTRK", "KEBAD", "CHI2", "BOTH", "DTH", "RATIO", "GOOD"};
   const int catColor[NCAT] = {kGray + 1, kGray + 2, kOrange + 7, kMagenta + 1, kRed + 1, kGreen + 2, kAzure + 2};
   const char *cfName[5] = {"hNoDth_", "hNoRat_", "hNoMatch_", "hNoChi2_", "hAny_"};
   const char *cfLbl[5] = {"noDth", "noRatio", "noMatch", "noChi2", "anyFit"};

   TH1D *gen[2] = {nullptr, nullptr}, *cat[2][NCAT] = {}, *cf[2][5] = {};
   for (int f = 0; f < 2; ++f)
      for (int s = s0; s <= s1; ++s) {
         TFile *fi = TFile::Open(TString::Format("%s/acc_split_%s_%s_s%d.root", dir.Data(), tag.Data(), fitName[f], s));
         if (!fi || fi->IsZombie()) {
            printf("\033[1;31mmissing %s seed %d\033[0m\n", fitName[f], s);
            continue;
         }
         auto add = [&](TH1D *&dst, const char *nm, const char *as) {
            auto *h = (TH1D *)fi->Get(nm);
            if (!h)
               return;
            if (!dst) {
               dst = (TH1D *)h->Clone(as);
               dst->SetDirectory(nullptr);
            } else
               dst->Add(h);
         };
         add(gen[f], "hGen_" + tag, TString::Format("gen_%s", fitName[f]));
         for (int c = 0; c < NCAT; ++c)
            add(cat[f][c], TString::Format("hCat_%s_%s", tag.Data(), catName[c]),
                TString::Format("cat_%s_%s", fitName[f], catName[c]));
         for (int k = 0; k < 5; ++k)
            add(cf[f][k], TString::Format("%s%s", cfName[k], tag.Data()),
                TString::Format("cf_%s_%d", fitName[f], k));
         fi->Close();
      }
   if (!gen[0] || !gen[1]) {
      printf("\033[1;31mnothing merged\033[0m\n");
      return;
   }

   TH1D *acc[2], *accCf[2][5];
   for (int f = 0; f < 2; ++f) {
      acc[f] = (TH1D *)cat[f][CAT_GOOD]->Clone(TString::Format("acc_%s", fitName[f]));
      acc[f]->Divide(cat[f][CAT_GOOD], gen[f], 1, 1, "B");
      for (int k = 0; k < 5; ++k) {
         accCf[f][k] = (TH1D *)cf[f][k]->Clone(TString::Format("accCf_%s_%d", fitName[f], k));
         accCf[f][k]->Divide(cf[f][k], gen[f], 1, 1, "B");
      }
   }

   for (int f = 0; f < 2; ++f) {
      printf("\n===== %s, %d seeds merged, %.0f generated reactions =====\n", fitLbl[f], s1 - s0 + 1,
             gen[f]->Integral());
      printf("  theta_cm     gen |");
      for (int c = NCAT - 1; c >= 0; --c)
         printf(" %6s", catName[c]);
      printf(" |  nominal");
      for (int k = 0; k < 5; ++k)
         printf(" %7s", cfLbl[k]);
      printf("\n");
      for (int b = 1; b <= gen[f]->GetNbinsX(); ++b) {
         double g = gen[f]->GetBinContent(b), ctr = gen[f]->GetBinCenter(b);
         if (g < 5 || ctr < 15 || ctr > 155)
            continue;
         printf("  %3.0f-%3.0f %6.0f |", gen[f]->GetBinLowEdge(b), gen[f]->GetBinLowEdge(b) + gen[f]->GetBinWidth(b),
                g);
         for (int c = NCAT - 1; c >= 0; --c)
            printf(" %5.1f%%", 100.0 * cat[f][c]->GetBinContent(b) / g);
         printf(" |   %6.3f", acc[f]->GetBinContent(b));
         for (int k = 0; k < 5; ++k)
            printf("  %6.3f", accCf[f][k]->GetBinContent(b));
         printf("\n");
      }
   }

   printf("\n===== the two fitters, side by side =====\n");
   printf("  theta_cm |   UKF acc  chi2-lost  no-chi2 |  GF acc  chi2-lost  no-chi2 |  GF/UKF nominal  GF/UKF no-chi2\n");
   for (int b = 1; b <= gen[0]->GetNbinsX(); ++b) {
      double ctr = gen[0]->GetBinCenter(b);
      if (ctr < 15 || ctr > 155 || gen[0]->GetBinContent(b) < 5)
         continue;
      double au = acc[0]->GetBinContent(b), ag = acc[1]->GetBinContent(b);
      double nu = accCf[0][3]->GetBinContent(b), ng = accCf[1][3]->GetBinContent(b);
      printf("  %3.0f-%3.0f  |   %6.3f    %5.1f%%    %6.3f |  %6.3f    %5.1f%%    %6.3f |     %6.3f          %6.3f\n",
             gen[0]->GetBinLowEdge(b), gen[0]->GetBinLowEdge(b) + gen[0]->GetBinWidth(b), au,
             100.0 * cat[0][CAT_CHI2]->GetBinContent(b) / gen[0]->GetBinContent(b), nu, ag,
             100.0 * cat[1][CAT_CHI2]->GetBinContent(b) / gen[1]->GetBinContent(b), ng,
             au > 0 ? ag / au : 0, nu > 0 ? ng / nu : 0);
   }

   TCanvas *c1 = new TCanvas("c1", "split, merged", 1500, 950);
   c1->Divide(2, 2);
   for (int f = 0; f < 2; ++f) {
      c1->cd(1 + f);
      auto *st = new THStack(TString::Format("st%d", f),
                             TString::Format("%s: where each reaction is lost;#theta_{cm} [deg];fraction", fitLbl[f]));
      for (int c = 0; c < NCAT; ++c) {
         auto *h = (TH1D *)cat[f][c]->Clone(TString::Format("fr_%s_%s", fitName[f], catName[c]));
         h->Divide(h, gen[f], 1, 1, "B");
         h->SetFillColor(catColor[c]);
         h->SetLineColor(catColor[c]);
         st->Add(h);
      }
      st->Draw("hist");
      st->GetXaxis()->SetRangeUser(15, 155);
      st->SetMaximum(1.05);
      auto *lg = new TLegend(0.60, 0.13, 0.89, 0.45);
      for (int c = NCAT - 1; c >= 0; --c) {
         auto *d = new TH1D(TString::Format("dm_%d_%s", f, catName[c]), "", 1, 0, 1);
         d->SetFillColor(catColor[c]);
         lg->AddEntry(d, catName[c], "f");
      }
      lg->Draw();
   }
   auto style = [](TH1D *h, int col, int mk, int ls = 1) {
      h->SetMarkerStyle(mk);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
      h->SetLineStyle(ls);
   };
   c1->cd(3);
   style(acc[0], kAzure + 2, 20);
   style(acc[1], kRed + 1, 21);
   style(accCf[0][3], kAzure + 2, 24, 2);
   style(accCf[1][3], kRed + 1, 25, 2);
   acc[0]->SetTitle("acceptance: nominal (solid) vs no #chi^{2} cut (open);#theta_{cm} [deg];acceptance");
   acc[0]->GetXaxis()->SetRangeUser(15, 155);
   acc[0]->GetYaxis()->SetRangeUser(0, 1.15);
   acc[0]->Draw("E1");
   acc[1]->Draw("E1 same");
   accCf[0][3]->Draw("E1 same");
   accCf[1][3]->Draw("E1 same");
   auto *lg3 = new TLegend(0.50, 0.13, 0.89, 0.38);
   lg3->AddEntry(acc[0], "UKF, all cuts", "lp");
   lg3->AddEntry(acc[1], "GENFIT, all cuts", "lp");
   lg3->AddEntry(accCf[0][3], "UKF, no #chi^{2} cut", "lp");
   lg3->AddEntry(accCf[1][3], "GENFIT, no #chi^{2} cut", "lp");
   lg3->Draw();

   c1->cd(4);
   auto *fu = (TH1D *)cat[0][CAT_CHI2]->Clone("fchi2_ukf");
   fu->Divide(cat[0][CAT_CHI2], gen[0], 1, 1, "B");
   auto *fg = (TH1D *)cat[1][CAT_CHI2]->Clone("fchi2_gf");
   fg->Divide(cat[1][CAT_CHI2], gen[1], 1, 1, "B");
   style(fu, kAzure + 2, 20);
   style(fg, kRed + 1, 21);
   fu->SetTitle("fraction lost to the #chi^{2}/ndf cut alone;#theta_{cm} [deg];fraction of generated");
   fu->GetXaxis()->SetRangeUser(15, 155);
   fu->GetYaxis()->SetRangeUser(0, 0.8);
   fu->Draw("E1");
   fg->Draw("E1 same");
   auto *lg4 = new TLegend(0.60, 0.72, 0.89, 0.87);
   lg4->AddEntry(fu, "UKF", "lp");
   lg4->AddEntry(fg, "GENFIT", "lp");
   lg4->Draw();

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString png = here + "/" + dir + "/acc_split_merged_" + tag + ".png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
