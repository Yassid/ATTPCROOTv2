/// @file exc_panels_C14.C
/// @brief The six excited-state angular distributions, one per panel.
///
/// Same numbers as exc_angdist_C14.C and exc8_angdist_C14.C produce; this only draws them
/// separately instead of overlaid, so each level can be read without the strong ones dominating
/// the frame.
///
/// THE AXIS MATTERS MORE THAN IT SOUNDS. With logy + common (the default) every panel shares one
/// logarithmic range, which is what the overlaid thesis figure uses -- relative strengths between
/// levels are then readable and the panels can be compared with each other and with that figure.
/// With logy = kFALSE and common = kFALSE each panel is stretched to its own linear range, which
/// makes a single shape easiest to read but makes panels look wildly different from one another
/// for reasons that are pure axis choice. Both are useful; they are not the same picture.
///
///   root -b -q 'exc_panels_C14.C()'

void exc_panels_C14(TString mFile = "plots/exc_angdist_gfex.root", TString uFile = "plots/exc8_angdist_hi.root",
                    Bool_t corrected = kTRUE, Bool_t logy = kTRUE, Bool_t common = kTRUE, TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const char *pre = corrected ? "dsdo" : "yield";

   struct P {
      const char *name, *file;
      int idx, col;
   };
   const std::vector<P> LV = {{"6.09 MeV  (1^{-})", nullptr, 0, kAzure + 2},
                              {"6.70 MeV  (3^{-})", nullptr, 1, kRed + 1},
                              {"7.00 MeV  (2^{+})", nullptr, 2, kGreen + 3},
                              {"7.27 MeV  (2^{-})", nullptr, 3, kOrange + 7},
                              {"8.53 MeV", nullptr, 0, kAzure + 2},
                              {"9.36 MeV  (blend)", nullptr, 1, kRed + 1}};

   TFile *fm = TFile::Open(here + "/" + mFile);
   TFile *fu = TFile::Open(here + "/" + uFile);
   if (!fm || fm->IsZombie() || !fu || fu->IsZombie()) {
      printf("\033[1;31mcannot open %s or %s\033[0m\n", mFile.Data(), uFile.Data());
      return;
   }

   TCanvas *c = new TCanvas("cp", "excited states", 1500, 900);
   c->Divide(3, 2);
   for (size_t i = 0; i < LV.size(); ++i) {
      TFile *f = (i < 4) ? fm : fu;
      auto *h = (TH1D *)f->Get(TString::Format("%s_%d", pre, LV[i].idx));
      if (!h) {
         printf("  %s: no %s_%d\n", LV[i].name, pre, LV[i].idx);
         continue;
      }
      auto *hc = (TH1D *)h->Clone(TString::Format("p%zu", i));
      hc->SetDirectory(nullptr);
      c->cd(i + 1);
      gPad->SetGridy();
      if (logy)
         gPad->SetLogy();
      hc->SetTitle(TString::Format("%s;#theta_{cm} [deg];%s", LV[i].name,
                                   corrected ? "d#sigma/d#Omega  [arb]" : "counts / bin"));
      hc->SetMarkerStyle(20);
      hc->SetMarkerSize(1.3);
      hc->SetMarkerColor(LV[i].col);
      hc->SetLineColor(LV[i].col);
      hc->SetLineWidth(2);
      double mx = 0;
      for (int b = 1; b <= hc->GetNbinsX(); ++b)
         mx = std::max(mx, hc->GetBinContent(b) + hc->GetBinError(b));
      if (common) {
         hc->SetMinimum(logy ? 3.0 : 0.0);
         hc->SetMaximum(logy ? 600.0 : 260.0);
      } else {
         hc->SetMinimum(logy ? std::max(1.0, mx * 0.02) : 0.0);
         hc->SetMaximum(mx * (logy ? 2.0 : 1.25));
      }
      hc->GetXaxis()->SetRangeUser(20, 140);
      hc->Draw("E1");

      // a bin that is empty because its fit was dropped is NOT a measurement of zero, and on a
      // linear axis it is indistinguishable from one. Mark it.
      for (int b = 1; b <= hc->GetNbinsX(); ++b) {
         if (hc->GetBinContent(b) > 0 || hc->GetBinCenter(b) < 20 || hc->GetBinCenter(b) > 140)
            continue;
         auto *tx = new TLatex(hc->GetBinCenter(b), common ? (logy ? 4.5 : mx * 0.06) : mx * 0.06, "#times");
         tx->SetTextColor(kGray + 2);
         tx->SetTextAlign(22);
         tx->SetTextSize(0.05);
         tx->Draw();
      }
   }
   TString png = here + "/plots/exc_panels_C14" + tag + ".png";
   c->SaveAs(png);
   printf("\n  (grey x = bin with no stored yield, i.e. dropped fit -- not a zero measurement)\n");
   printf("wrote %s\n\n", png.Data());
}
