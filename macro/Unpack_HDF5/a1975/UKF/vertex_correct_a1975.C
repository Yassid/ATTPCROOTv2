/// @file vertex_correct_a1975.C
/// @brief Vertex-z correction of the 16C Ex spectrum (resolution improvement).
///
/// The beam loses energy in the gas before the vertex, so the apparent Ex of the
/// elastic peak drifts with vertex-z. This fits the ELASTIC peak position per
/// vertex-z slice (mu(z)) and aligns every slice to Ex=0 (Ex_corr = Ex - mu(z)),
/// removing the vertex-z broadening. Reads the proton_kin.root cache (fast, no FRIB)
/// and compares the elastic-peak FWHM before vs after.
///
///   root -b -q 'vertex_correct_a1975.C()'

#include <vector>

void vertex_correct_a1975(TString cache = "proton_kin.root", int nz = 28, double zlo = 60, double zhi = 1000,
                         double elasticWin = 4.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *f = TFile::Open(cache);
   TNtuple *pk = (TNtuple *)f->Get("pk");
   float ke, theta, vertexz, thcm, ex, chi2ndf;
   pk->SetBranchAddress("ke", &ke);
   pk->SetBranchAddress("vertexz", &vertexz);
   pk->SetBranchAddress("ex", &ex);
   Long64_t n = pk->GetEntries();

   // --- pass 1: fit elastic-peak position per vertex-z slice ---
   double dz = (zhi - zlo) / nz;
   std::vector<TH1F *> slice(nz);
   for (int s = 0; s < nz; ++s)
      slice[s] = new TH1F(Form("sl%d", s), "", 120, -elasticWin, elasticWin);
   for (Long64_t i = 0; i < n; ++i) {
      pk->GetEntry(i);
      if (std::fabs(ex) > elasticWin)
         continue;
      int s = (int)((vertexz - zlo) / dz);
      if (s >= 0 && s < nz)
         slice[s]->Fill(ex);
   }
   TGraph *gmu = new TGraph();
   int np = 0;
   for (int s = 0; s < nz; ++s) {
      if (slice[s]->GetEntries() < 40)
         continue;
      double pk0 = slice[s]->GetBinCenter(slice[s]->GetMaximumBin());
      TF1 g("g", "gaus", pk0 - 1.5, pk0 + 1.5);
      g.SetParameters(slice[s]->GetMaximum(), pk0, 0.6);
      slice[s]->Fit(&g, "QRN");
      gmu->SetPoint(np++, zlo + (s + 0.5) * dz, g.GetParameter(1));
   }
   // smooth correction mu(z) = pol2 fit
   TF1 *fmu = new TF1("fmu", "pol2", zlo, zhi);
   gmu->Fit(fmu, "QRN");

   // --- pass 2: apply Ex_corr = Ex - mu(z); compare elastic peaks ---
   TH1F *hraw = new TH1F("hraw", "16C E_{x}: raw vs vertex-corrected;E_{x} [MeV];protons", 240, -5, 10);
   TH1F *hcor = new TH1F("hcor", "", 240, -5, 10);
   for (Long64_t i = 0; i < n; ++i) {
      pk->GetEntry(i);
      double mu = fmu->Eval(std::min(std::max((double)vertexz, zlo), zhi));
      hraw->Fill(ex);
      hcor->Fill(ex - mu);
   }
   auto fwhm = [](TH1F *h) {
      TF1 g("g", "gaus", -2, 2);
      g.SetParameters(h->GetMaximum(), 0, 0.8);
      h->Fit(&g, "QRN");
      return std::make_pair(2.3548 * g.GetParameter(2), g.GetParameter(1));
   };
   auto [fw_raw, mu_raw] = fwhm(hraw);
   auto [fw_cor, mu_cor] = fwhm(hcor);
   printf("\nElastic peak  RAW: mu=%.3f FWHM=%.3f MeV   CORRECTED: mu=%.3f FWHM=%.3f MeV  (%.0f%% narrower)\n", mu_raw,
          fw_raw, mu_cor, fw_cor, 100.0 * (fw_raw - fw_cor) / fw_raw);
   printf("mu(z) correction = %.4g + %.4g*z + %.4g*z^2  (MeV, z in mm)\n", fmu->GetParameter(0), fmu->GetParameter(1),
          fmu->GetParameter(2));

   TCanvas *c = new TCanvas("c", "vcorr", 1300, 560);
   c->Divide(2, 1);
   c->cd(1);
   gmu->SetTitle("elastic peak position vs vertex z;vertex z [mm];E_{x} peak [MeV]");
   gmu->SetMarkerStyle(20);
   gmu->Draw("AP");
   fmu->SetLineColor(kRed);
   fmu->Draw("same");
   c->cd(2);
   hraw->SetLineColor(kGray + 2);
   hraw->SetLineWidth(2);
   hcor->SetLineColor(kRed + 1);
   hcor->SetLineWidth(2);
   hraw->Draw("hist");
   hcor->Draw("hist same");
   TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg->AddEntry(hraw, Form("raw (FWHM %.2f)", fw_raw), "l");
   leg->AddEntry(hcor, Form("corrected (FWHM %.2f)", fw_cor), "l");
   leg->Draw();
   c->SaveAs("vertex_correction.png");
   printf("saved vertex_correction.png\n");
}
