/// @file ex_slices_C14.C
/// @brief Ex spectrum in theta_cm slices, UKF vs GENFIT -- is the elastic peak still there?
///
/// FRESCO says the elastic cross section at theta_cm 70-110 deg is only ~4x below its value at
/// 40-45 deg, yet both fitters deliver ~10x less yield than predicted in an Ex window centred on
/// zero. Two possibilities: the events are gone, or they are present at the wrong Ex. The
/// per-bin totals over ALL Ex are equal between fitters and do not collapse, which points at the
/// second. This macro looks directly: one Ex spectrum per theta_cm slice, and a gaussian fit to
/// whatever peak is there, so the elastic locus can be tracked rather than assumed at zero.
///
/// The kinematic origin of a drift: at these theta_cm the recoil proton comes out at small
/// theta_lab (theta_lab ~ (180-theta_cm)/2, so 35-55 deg) with 11-30 MeV. dEx/dtheta_lab is
/// largest exactly there, so a small angular bias shows up as a large Ex shift -- which is why
/// this has to be checked before any conclusion about either fitter.
///
///   root -b -q 'ex_slices_C14.C()'

void ex_slices_C14(Double_t exLo = -3.0, Double_t exHi = 3.0, Int_t nb = 120)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const int NS = 2;
   const char *file[NS] = {"plots/proton_kin_300_ukf_nc.root", "plots/proton_kin_300gfx_nc.root"};
   const char *lbl[NS] = {"UKF", "GENFIT"};
   const int col[NS] = {kAzure + 2, kRed + 1};

   const int NSL = 12;
   const double slo[NSL] = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130};
   const double shi[NSL] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 145};

   TTree *t[NS];
   TFile *f[NS];
   for (int i = 0; i < NS; ++i) {
      f[i] = TFile::Open(here + "/" + file[i]);
      if (!f[i] || f[i]->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", file[i]);
         return;
      }
      t[i] = (TTree *)f[i]->Get("pk");
   }

   TCanvas *c1 = new TCanvas("c1", "Ex slices", 1600, 1100);
   c1->Divide(4, 3);
   printf("\n===== elastic locus per theta_cm slice (gaussian on the tallest structure) =====\n");
   printf("  theta_cm |        UKF: mu    sigma      N |     GENFIT: mu    sigma      N\n");
   double mu[NS][NSL], sg[NS][NSL];
   for (int s = 0; s < NSL; ++s) {
      c1->cd(s + 1);
      double ymax = 0;
      TH1D *h[NS];
      for (int i = 0; i < NS; ++i) {
         h[i] = new TH1D(TString::Format("h%d_%d", i, s), "", nb, exLo, exHi);
         t[i]->Draw(TString::Format("ex>>h%d_%d", i, s), TString::Format("thcm>=%g&&thcm<%g", slo[s], shi[s]), "goff");
         h[i]->SetDirectory(nullptr);
         h[i]->SetLineColor(col[i]);
         h[i]->SetLineWidth(2);
         ymax = std::max(ymax, h[i]->GetMaximum());
         // fit the tallest structure in the slice
         int bm = 0;
         double vm = -1;
         for (int b = 1; b <= nb; ++b)
            if (h[i]->GetBinContent(b) > vm) {
               vm = h[i]->GetBinContent(b);
               bm = b;
            }
         double c0 = h[i]->GetBinCenter(bm);
         TF1 g(TString::Format("g%d_%d", i, s), "gaus", c0 - 0.5, c0 + 0.5);
         if (h[i]->Integral() > 30) {
            h[i]->Fit(&g, "QNR");
            mu[i][s] = g.GetParameter(1);
            sg[i][s] = std::fabs(g.GetParameter(2));
         } else {
            mu[i][s] = c0;
            sg[i][s] = 0;
         }
      }
      h[0]->SetTitle(TString::Format("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", slo[s], shi[s]));
      h[0]->SetMaximum(ymax * 1.25);
      h[0]->Draw("hist");
      h[1]->Draw("hist same");
      auto *l0 = new TLine(0, 0, 0, ymax * 1.25);
      l0->SetLineStyle(2);
      l0->SetLineColor(kGray + 2);
      l0->Draw();
      if (s == 0) {
         auto *lg = new TLegend(0.60, 0.72, 0.92, 0.90);
         for (int i = 0; i < NS; ++i)
            lg->AddEntry(h[i], lbl[i], "l");
         lg->Draw();
      }
      printf("  %3.0f-%3.0f |  %+12.3f %8.3f %6.0f |  %+12.3f %8.3f %6.0f\n", slo[s], shi[s], mu[0][s], sg[0][s],
             h[0]->Integral(), mu[1][s], sg[1][s], h[1]->Integral());
   }

   TString png = here + "/plots/ex_slices_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n", png.Data());

   // the drift curve itself, for use as a ridge-tracking window
   TCanvas *c2 = new TCanvas("c2", "elastic locus", 900, 650);
   auto *gU = new TGraph(), *gG = new TGraph();
   for (int s = 0; s < NSL; ++s) {
      gU->SetPoint(s, 0.5 * (slo[s] + shi[s]), mu[0][s]);
      gG->SetPoint(s, 0.5 * (slo[s] + shi[s]), mu[1][s]);
   }
   gU->SetMarkerStyle(20);
   gU->SetMarkerColor(kAzure + 2);
   gU->SetLineColor(kAzure + 2);
   gU->SetLineWidth(2);
   gG->SetMarkerStyle(21);
   gG->SetMarkerColor(kRed + 1);
   gG->SetLineColor(kRed + 1);
   gG->SetLineWidth(2);
   gU->SetTitle("elastic locus vs #theta_{cm};#theta_{cm} [deg];fitted E_{x} peak [MeV]");
   gU->GetYaxis()->SetRangeUser(-2.5, 1.0);
   gU->Draw("ALP");
   gG->Draw("LP same");
   auto *z = new TLine(20, 0, 145, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   auto *lg2 = new TLegend(0.62, 0.75, 0.89, 0.89);
   lg2->AddEntry(gU, "UKF", "lp");
   lg2->AddEntry(gG, "GENFIT", "lp");
   lg2->Draw();
   TString png2 = here + "/plots/ex_locus_C14.png";
   c2->SaveAs(png2);
   printf("wrote %s\n\n", png2.Data());
}
