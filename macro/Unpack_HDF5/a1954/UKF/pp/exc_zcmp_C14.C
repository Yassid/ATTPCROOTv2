/// @file exc_zcmp_C14.C
/// @brief The six excited-state angular distributions with and without the vertex-z window.
///
/// The window exists so that the excited-state yields and the elastic normalisation are counted
/// over the same target thickness. It costs statistics, so the question this plot answers is
/// whether it also costs, or changes, the SHAPE -- which is what the multipole comparison uses.
///
/// The two are drawn on a common scale, and each is separately normalised to its own integral in
/// the lower row so that shape can be compared independently of the rate. If the window only
/// removed events, the normalised curves lie on top of each other; where they separate, the
/// window is selecting a different mixture rather than a smaller sample of the same one.
///
///   root -b -q 'exc_zcmp_C14.C()'

void exc_zcmp_C14(TString mZ = "plots/exc_angdist_gfex.root", TString uZ = "plots/exc8_angdist_hi.root",
                  TString mN = "plots/exc_angdist_noz.root", TString uN = "plots/exc8_angdist_noz.root",
                  TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   struct P {
      const char *name;
      int idx, col;
      bool upper;
   };
   const std::vector<P> LV = {{"6.09 MeV  (1^{-})", 0, kAzure + 2, false}, {"6.70 MeV  (3^{-})", 1, kRed + 1, false},
                              {"7.00 MeV  (2^{+})", 2, kGreen + 3, false}, {"7.27 MeV  (2^{-})", 3, kOrange + 7, false},
                              {"8.53 MeV", 0, kAzure + 2, true},          {"9.36 MeV  (blend)", 1, kRed + 1, true}};

   TFile *fmZ = TFile::Open(here + "/" + mZ), *fuZ = TFile::Open(here + "/" + uZ);
   TFile *fmN = TFile::Open(here + "/" + mN), *fuN = TFile::Open(here + "/" + uN);
   if (!fmZ || fmZ->IsZombie() || !fuZ || fuZ->IsZombie() || !fmN || fmN->IsZombie() || !fuN || fuN->IsZombie()) {
      printf("\033[1;31mmissing one of the four inputs\033[0m\n");
      return;
   }

   TCanvas *c = new TCanvas("cz", "z window comparison", 1550, 950);
   c->Divide(3, 2);
   printf("\n  ratio (z window / no window) of the integral, and of the shape at the extremes\n");
   printf("  level        | integral | fwd bin | bwd bin\n");

   for (size_t i = 0; i < LV.size(); ++i) {
      auto *hZ = (TH1D *)(LV[i].upper ? fuZ : fmZ)->Get(TString::Format("dsdo_%d", LV[i].idx));
      auto *hN = (TH1D *)(LV[i].upper ? fuN : fmN)->Get(TString::Format("dsdo_%d", LV[i].idx));
      if (!hZ || !hN) {
         printf("  %-12s missing\n", LV[i].name);
         continue;
      }
      auto *a = (TH1D *)hZ->Clone(TString::Format("z%zu", i));
      auto *b = (TH1D *)hN->Clone(TString::Format("n%zu", i));
      a->SetDirectory(nullptr);
      b->SetDirectory(nullptr);

      c->cd(i + 1);
      gPad->SetLogy();
      gPad->SetGridy();
      b->SetTitle(TString::Format("%s;#theta_{cm} [deg];d#sigma/d#Omega  [arb]", LV[i].name));
      b->SetMinimum(3);
      b->SetMaximum(600);
      b->GetXaxis()->SetRangeUser(20, 140);
      // no window: open grey squares, drawn first so the window points sit on top
      b->SetMarkerStyle(25);
      b->SetMarkerSize(1.3);
      b->SetMarkerColor(kGray + 2);
      b->SetLineColor(kGray + 2);
      b->SetLineWidth(2);
      b->Draw("E1");
      a->SetMarkerStyle(20);
      a->SetMarkerSize(1.3);
      a->SetMarkerColor(LV[i].col);
      a->SetLineColor(LV[i].col);
      a->SetLineWidth(2);
      a->Draw("E1 same");

      if (i == 0) {
         auto *lg = new TLegend(0.15, 0.15, 0.62, 0.32);
         lg->AddEntry(b, "no z window", "lp");
         lg->AddEntry(a, "z 10-400 mm", "lp");
         lg->SetTextSize(0.045);
         lg->Draw();
      }

      // how much of the change is rate and how much is shape
      double iz = a->Integral(), in = b->Integral();
      int fwd = -1, bwd = -1;
      for (int k = 1; k <= a->GetNbinsX(); ++k) {
         if (a->GetBinContent(k) > 0 && b->GetBinContent(k) > 0) {
            if (fwd < 0)
               fwd = k;
            bwd = k;
         }
      }
      auto r = [&](int k) { return (k > 0 && b->GetBinContent(k) > 0) ? a->GetBinContent(k) / b->GetBinContent(k) : 0; };
      printf("  %-12s |  %6.3f  |  %6.3f |  %6.3f\n", LV[i].name, in > 0 ? iz / in : 0, r(fwd), r(bwd));
   }

   TString png = here + "/plots/exc_zcmp_C14" + tag + ".png";
   c->SaveAs(png);
   printf("\n  A ratio that is the same in the forward and backward bins means the window only\n"
          "  removed events; one that differs means it changed the angular mixture.\n");
   printf("wrote %s\n\n", png.Data());
}
