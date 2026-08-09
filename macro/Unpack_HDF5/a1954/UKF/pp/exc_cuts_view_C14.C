/// @file exc_cuts_view_C14.C
/// @brief Draw what the excited-state fit is actually selecting, on the Ex vs theta_cm plane.
///
/// The per-angle fit (exc_angdist_C14.C) holds the four component centroids FIXED at the
/// calibrated level energies. That is only legitimate if the measured ridges do not move with
/// angle. This macro puts the fixed centroids, the fit range and the angular binning on top of
/// the data so the assumption can be checked rather than trusted, and overlays the ridge position
/// measured bin by bin.
///
///   root -b -q 'exc_cuts_view_C14.C()'

void exc_cuts_view_C14(TString cache = "plots/proton_kin_300gfx_ex.root", Double_t gain = 1.078,
                       TString fitFile = "plots/exc_angdist_gfex.root",
                       Double_t cmMin = 20.0, Double_t cmMax = 140.0, Double_t dcm = 10.0,
                       Double_t fitLo = 5.4, Double_t fitHi = 8.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const int NLV = 4;
   const double MU0[NLV] = {6.09, 6.70, 7.00, 7.27};
   const char *LVNAME[NLV] = {"6.09 (1^{-})", "6.70 (3^{-})", "7.00 (2^{+})", "7.27 (2^{-})"};
   const int LVCOL[NLV] = {kAzure + 2, kRed + 1, kGreen + 3, kOrange + 7};
   double mu[NLV];
   for (int i = 0; i < NLV; ++i)
      mu[i] = 6.094 + gain * (MU0[i] - 6.094);

   // the per-bin shift the fit actually applied
   TGraph *gsh = nullptr;
   { TFile *ff = TFile::Open(here + "/" + fitFile);
     if (ff && !ff->IsZombie()) { auto *g = (TGraph *)ff->Get("shift");
       if (g) gsh = (TGraph *)g->Clone("gsh"); ff->Close(); } }
   if (!gsh)
      printf("\033[1;33mno shift found in %s -- drawing the UNSHIFTED centroids only\033[0m\n", fitFile.Data());

   TFile *f = TFile::Open(here + "/" + cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");

   TCanvas *c = new TCanvas("c", "what the fit selects", 1500, 620);
   c->Divide(2, 1);

   // ---- left: the plane with the fixed centroids and the fit window
   c->cd(1);
   gPad->SetLogz();
   auto *h2 = new TH2D("h2", "E_{x} vs #theta_{cm}: fixed fit centroids and the fit range;"
                             "#theta_{cm} [deg];E_{x} [MeV]",
                       60, 0, 180, 160, 4.5, 9.0);
   t->Draw("ex:thcm>>h2", "", "goff");
   h2->SetDirectory(nullptr);
   h2->Draw("colz");
   for (int i = 0; i < NLV; ++i) {
      auto *l = new TLine(cmMin, mu[i], cmMax, mu[i]);
      l->SetLineColor(LVCOL[i]);
      l->SetLineWidth(2);
      l->Draw();
      auto *tx = new TLatex(cmMax + 2, mu[i] - 0.05, LVNAME[i]);
      tx->SetTextColor(LVCOL[i]);
      tx->SetTextSize(0.028);
      tx->Draw();
   }
   for (double y : {fitLo, fitHi}) { // the fit range
      auto *l = new TLine(cmMin, y, cmMax, y);
      l->SetLineColor(kWhite);
      l->SetLineStyle(2);
      l->SetLineWidth(2);
      l->Draw();
   }
   for (double x = cmMin; x <= cmMax + 0.1; x += dcm) { // the angular binning
      auto *l = new TLine(x, fitLo, x, fitHi);
      l->SetLineColor(kGray + 1);
      l->SetLineStyle(3);
      l->Draw();
   }

   // ---- right: the same, zoomed, with the MEASURED ridge per bin on top
   c->cd(2);
   gPad->SetLogz();
   auto *h3 = new TH2D("h3", "zoom: centroids AS FITTED (solid) vs unshifted (dotted) vs measured ridge (white);"
                             "#theta_{cm} [deg];E_{x} [MeV]",
                       (int)((cmMax - cmMin) / 2.5), cmMin, cmMax, 120, 5.4, 8.2);
   t->Draw("ex:thcm>>h3", "", "goff");
   h3->SetDirectory(nullptr);
   h3->Draw("colz");
   // the centroids AS FITTED: flat position plus the per-bin shift
   for (int i = 0; i < NLV; ++i) {
      auto *gc = new TGraph();
      int m = 0;
      for (double lo = cmMin; lo < cmMax - 0.1; lo += dcm) {
         double x = lo + dcm / 2;
         double sh = gsh ? gsh->Eval(x) : 0.0;
         gc->SetPoint(m++, x, mu[i] + sh);
      }
      gc->SetLineColor(LVCOL[i]);
      gc->SetLineWidth(3);
      gc->Draw("L same");
      auto *fl = new TLine(cmMin, mu[i], cmMax, mu[i]); // unshifted, for reference
      fl->SetLineColor(LVCOL[i]);
      fl->SetLineStyle(3);
      fl->Draw();
   }
   // measured ridge: mode of the Ex spectrum in each angular bin, restricted to the fit range
   auto *gr = new TGraph();
   int n = 0;
   printf("\n  theta_cm |  measured ridge  |  nearest fixed centroid  |  offset\n");
   for (double lo = cmMin; lo < cmMax - 0.1; lo += dcm) {
      auto *hx = new TH1D(TString::Format("hr%d", (int)lo), "", 56, fitLo, fitHi);
      t->Draw(TString::Format("ex>>hr%d", (int)lo), TString::Format("thcm>=%g&&thcm<%g", lo, lo + dcm), "goff");
      hx->SetDirectory(nullptr);
      if (hx->Integral() > 60) {
         hx->Smooth(2);
         double pk = hx->GetBinCenter(hx->GetMaximumBin());
         gr->SetPoint(n++, lo + dcm / 2, pk);
         int best = 0;
         for (int i = 1; i < NLV; ++i)
            if (std::fabs(pk - mu[i]) < std::fabs(pk - mu[best]))
               best = i;
         printf("  %3.0f-%3.0f  |      %6.3f      |   %-14s %6.3f |  %+6.3f\n", lo, lo + dcm, pk, LVNAME[best],
                mu[best], pk - mu[best]);
      }
      delete hx;
   }
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kWhite);
   gr->SetMarkerSize(1.3);
   gr->SetLineColor(kWhite);
   gr->SetLineWidth(2);
   gr->Draw("LP same");

   TString png = here + "/plots/exc_cuts_view_C14.png";
   c->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
