/// @file exsel_view_C14.C
/// @brief What the excited-state analysis actually selects, drawn on the Ex vs theta_cm plane.
///
/// Both extractions hold their component centroids FIXED at the calibrated level energies and fit
/// only amplitudes. That is legitimate only if the measured ridges sit where the fit puts them and
/// do not move with angle, so this draws the fixed positions, the fit ranges and the angular
/// binning on top of the data, with the ridge measured bin by bin for comparison.
///
/// The vertex-z window is applied, because it is part of the selection: these are the events the
/// yields are extracted from, not the whole data set.
///
/// Panels: the plane over the full range, then the two fitted regions zoomed -- the 6-7.5 MeV
/// multiplet (fitted 5.4-8.0) and the 8.5/9.4 structures (fitted 7.9-10.3).
///
///   root -b -q 'exsel_view_C14.C()'

void exsel_view_C14(TString cache = "plots/proton_kin_300gfx_ex.root", Double_t zMin = 10.0, Double_t zMax = 400.0,
                    Double_t gain = 1.078, Double_t cmMin = 20.0, Double_t cmMax = 140.0, TString tag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *f = TFile::Open(here + "/" + cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t)
      return;
   TString zc = TString::Format("vertexz>%g&&vertexz<%g", zMin, zMax);

   // the multiplet centroids are the literature energies mapped through the calibration gain;
   // the two upper ones are as measured on the angle-integrated spectrum
   const int NM = 4;
   const double MU0[NM] = {6.09, 6.70, 7.00, 7.27};
   const char *MNAME[NM] = {"6.09 (1^{-})", "6.70 (3^{-})", "7.00 (2^{+})", "7.27 (2^{-})"};
   const int MCOL[NM] = {kAzure + 2, kRed + 1, kGreen + 3, kOrange + 7};
   double mu[NM];
   for (int i = 0; i < NM; ++i)
      mu[i] = 6.094 + gain * (MU0[i] - 6.094);
   const int NU = 2;
   const double MUU[NU] = {8.533, 9.363};
   const char *UNAME[NU] = {"8.53", "9.36"};
   const int UCOL[NU] = {kAzure + 2, kRed + 1};

   TCanvas *c = new TCanvas("cs", "selection", 1550, 950);
   c->Divide(2, 2);

   auto drawPlane = [&](int pad, double eLo, double eHi, int nEx, const char *title) {
      c->cd(pad);
      gPad->SetLogz();
      auto *h = new TH2D(TString::Format("hp%d", pad), TString::Format("%s;#theta_{cm} [deg];E_{x} [MeV]", title),
                         (int)((cmMax - cmMin) / 2.5), cmMin, cmMax, nEx, eLo, eHi);
      t->Draw(TString::Format("ex:thcm>>hp%d", pad), zc, "goff");
      h->SetDirectory(nullptr);
      h->Draw("colz");
      return h;
   };

   // ---- 1. the whole plane, both fitted regions marked
   drawPlane(1, 4.5, 10.8, 130, "E_{x} vs #theta_{cm}: both fit ranges");
   for (double y : {5.4, 8.0}) { // multiplet fit range
      auto *l = new TLine(cmMin, y, cmMax, y);
      l->SetLineColor(kWhite);
      l->SetLineWidth(2);
      l->SetLineStyle(2);
      l->Draw();
   }
   for (double y : {7.9, 10.3}) { // upper-states fit range
      auto *l = new TLine(cmMin, y, cmMax, y);
      l->SetLineColor(kMagenta);
      l->SetLineWidth(2);
      l->SetLineStyle(3);
      l->Draw();
   }

   // ---- 2. the multiplet, with the measured ridge
   drawPlane(2, 5.3, 8.1, 112, "multiplet: fixed centroids (lines) vs measured ridge (white)");
   for (int i = 0; i < NM; ++i) {
      auto *l = new TLine(cmMin, mu[i], cmMax, mu[i]);
      l->SetLineColor(MCOL[i]);
      l->SetLineWidth(2);
      l->Draw();
      auto *tx = new TLatex(cmMax + 1, mu[i] - 0.04, MNAME[i]);
      tx->SetTextColor(MCOL[i]);
      tx->SetTextSize(0.030);
      tx->Draw();
   }
   {
      auto *gr = new TGraph();
      int n = 0;
      printf("\n  MULTIPLET: measured ridge vs the nearest fixed centroid\n");
      printf("  theta_cm | ridge | nearest        | offset\n");
      for (double lo = cmMin; lo < cmMax - 0.1; lo += 10) {
         auto *hx = new TH1D(TString::Format("hm%d", (int)lo), "", 56, 5.4, 8.0);
         t->Draw(TString::Format("ex>>hm%d", (int)lo),
                 TString::Format("thcm>=%g&&thcm<%g&&", lo, lo + 10) + zc, "goff");
         hx->SetDirectory(nullptr);
         if (hx->Integral() > 40) {
            hx->Smooth(2);
            double pk = hx->GetBinCenter(hx->GetMaximumBin());
            gr->SetPoint(n++, lo + 5, pk);
            int best = 0;
            for (int i = 1; i < NM; ++i)
               if (std::fabs(pk - mu[i]) < std::fabs(pk - mu[best]))
                  best = i;
            printf("  %3.0f-%3.0f  | %5.3f | %-14s | %+6.3f\n", lo, lo + 10, pk, MNAME[best], pk - mu[best]);
         }
         delete hx;
      }
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kWhite);
      gr->SetLineColor(kWhite);
      gr->SetLineWidth(2);
      gr->SetMarkerSize(1.2);
      gr->Draw("LP same");
   }

   // ---- 3. the upper structures
   drawPlane(3, 7.7, 10.5, 112, "8.5 / 9.4 structures: fixed centroids");
   for (int i = 0; i < NU; ++i) {
      auto *l = new TLine(cmMin, MUU[i], cmMax, MUU[i]);
      l->SetLineColor(UCOL[i]);
      l->SetLineWidth(2);
      l->Draw();
      auto *tx = new TLatex(cmMax + 1, MUU[i] - 0.04, UNAME[i]);
      tx->SetTextColor(UCOL[i]);
      tx->SetTextSize(0.030);
      tx->Draw();
   }

   // ---- 4. the projections the fits actually see, in three angular slices
   c->cd(4);
   gPad->SetLogy();
   const int NS = 3;
   const double SLO[NS] = {20, 60, 100}, SHI[NS] = {60, 100, 140};
   const int SCOL[NS] = {kAzure + 2, kRed + 1, kGreen + 3};
   auto *lg = new TLegend(0.58, 0.70, 0.89, 0.88);
   double mx = 0;
   TH1D *hs[NS];
   for (int i = 0; i < NS; ++i) {
      hs[i] = new TH1D(TString::Format("hs%d", i), "E_{x} in angular slices;E_{x} [MeV];counts", 130, 4.5, 10.8);
      t->Draw(TString::Format("ex>>hs%d", i), TString::Format("thcm>=%g&&thcm<%g&&", SLO[i], SHI[i]) + zc, "goff");
      hs[i]->SetDirectory(nullptr);
      mx = std::max(mx, hs[i]->GetMaximum());
   }
   for (int i = 0; i < NS; ++i) {
      hs[i]->SetLineColor(SCOL[i]);
      hs[i]->SetLineWidth(2);
      hs[i]->SetMaximum(mx * 2.5);
      hs[i]->SetMinimum(0.5);
      hs[i]->Draw(i ? "hist same" : "hist");
      lg->AddEntry(hs[i], TString::Format("#theta_{cm} %.0f-%.0f", SLO[i], SHI[i]), "l");
   }
   for (int i = 0; i < NM; ++i) {
      auto *l = new TLine(mu[i], 0.5, mu[i], mx * 2.5);
      l->SetLineColor(MCOL[i]);
      l->SetLineStyle(2);
      l->Draw();
   }
   for (int i = 0; i < NU; ++i) {
      auto *l = new TLine(MUU[i], 0.5, MUU[i], mx * 2.5);
      l->SetLineColor(kGray + 2);
      l->SetLineStyle(2);
      l->Draw();
   }
   lg->Draw();

   TString png = here + "/plots/exsel_view_C14" + tag + ".png";
   c->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
