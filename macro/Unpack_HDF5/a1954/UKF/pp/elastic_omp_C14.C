/// @file elastic_omp_C14.C
/// @brief Which optical potential does the MEASURED elastic angular distribution prefer?
///
/// The extraction is elastic_wtrack_C14.C verbatim (width-tracking window mu +- k*w, no background
/// subtraction). What changes is the comparison: instead of one FRESCO/KD03 curve, the measured
/// shape is confronted with every global proton OMP in the PtolemyCpp library, each allowed ONE
/// free scale. That scale is the luminosity, which is unknown a priori, so only the SHAPE can
/// discriminate -- and the shape is what the dip position lives in.
///
/// Ptolemy and FRESCO were shown to agree to 0.05% on identical KD03 input, so the Ptolemy curves
/// are interchangeable with FRESCO ones here.
///
///   root -b -q 'elastic_omp_C14.C()'

namespace eo {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b)
      if (b > 0 && a >= 5 && a <= 175) g->SetPoint(g->GetN(), a, b);
   if (!g->GetN()) printf("\033[1;31mempty %s\033[0m\n", f.Data());
   return g;
}
} // namespace eo

void elastic_omp_C14(TString cache = "plots/proton_kin_cat5_s013.root",
                     TString accDir = "/mnt/f/a1954_C14_acc_catima_z10_490/", Double_t kSig = 2.5,
                     Double_t cmMin = 15, Double_t cmMax = 150, Double_t dcm = 5.0,
                     Double_t zMin = 10, Double_t zMax = 490, Double_t chi2Cut = 5.0,
                     Double_t fitLo = 20, Double_t fitHi = 148, TString tag = "omp")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = fd && !fd->IsZombie() ? (TNtuple *)fd->Get("pk") : nullptr;
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   TH1D *acc = fa && !fa->IsZombie() ? (TH1D *)fa->Get("hAcc_gs_sum") : nullptr;
   if (!t || !acc) { printf("\033[1;31mmissing cache or acceptance\033[0m\n"); return; }

   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   std::vector<double> mu(NB, 0), w(NB, 0), ctr(NB, 0);
   std::vector<bool> ok(NB, false);
   std::vector<TH1D *> hs(NB);
   float *v;
   for (int b = 0; b < NB; ++b) {
      hs[b] = new TH1D(Form("ho%d", b), "", 200, -6, 4);
      ctr[b] = cmMin + (b + 0.5) * dcm;
   }
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); v = t->GetArgs();
      if (v[5] > chi2Cut || v[0] <= 0 || v[2] < zMin || v[2] > zMax) continue;
      int b = (int)((v[3] - cmMin) / dcm);
      if (b >= 0 && b < NB) hs[b]->Fill(v[4]);
   }
   // ---- pass 1: locus and width (identical to elastic_wtrack_C14.C) --------------------------
   for (int b = 0; b < NB; ++b) {
      if (hs[b]->Integral() < 60) continue;
      TH1D *h = (TH1D *)hs[b]->Clone(Form("so%d", b)); h->Smooth(2);
      double pk = h->GetBinCenter(h->GetMaximumBin()), ymax = h->GetMaximum();
      int bm = h->GetMaximumBin(), lo = bm, hi = bm;
      while (lo > 1 && h->GetBinContent(lo) > 0.5 * ymax) --lo;
      while (hi < h->GetNbinsX() && h->GetBinContent(hi) > 0.5 * ymax) ++hi;
      double fwhm = h->GetBinCenter(hi) - h->GetBinCenter(lo);
      if (fwhm <= 0) continue;
      double ww = fwhm / 2.355, m = pk;
      for (int it = 0; it < 5; ++it) {
         double s = 0, n = 0;
         for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double x = h->GetBinCenter(i);
            if (std::fabs(x - m) > 1.2 * ww) continue;
            s += x * h->GetBinContent(i); n += h->GetBinContent(i);
         }
         if (n > 0) m = s / n;
      }
      mu[b] = m; w[b] = ww; ok[b] = true;
   }
   for (int b = 0; b < NB; ++b) if (!ok[b]) {
      int l = b, r = b;
      while (l >= 0 && !ok[l]) --l;
      while (r < NB && !ok[r]) ++r;
      if (l < 0 && r >= NB) continue;
      if (l < 0) { mu[b] = mu[r]; w[b] = w[r]; }
      else if (r >= NB) { mu[b] = mu[l]; w[b] = w[l]; }
      else { double f = (double)(b - l) / (r - l);
             mu[b] = mu[l] + f * (mu[r] - mu[l]); w[b] = w[l] + f * (w[r] - w[l]); }
   }
   std::vector<double> ms = mu, ws = w;
   for (int b = 1; b < NB - 1; ++b) { ms[b] = (mu[b-1]+mu[b]+mu[b+1])/3; ws[b] = (w[b-1]+w[b]+w[b+1])/3; }

   // ---- pass 2: count -----------------------------------------------------------------------
   std::vector<double> X, Y, E;
   printf("\n  extraction: window mu +- %.1f w, z %.0f-%.0f, chi2 < %.1f\n", kSig, zMin, zMax, chi2Cut);
   printf("\n  theta_cm   counts    acc     dsdo [arb]   err\n");
   for (int b = 0; b < NB; ++b) {
      double lo = ms[b] - kSig * ws[b], hi = ms[b] + kSig * ws[b];
      double y = hs[b]->Integral(hs[b]->FindBin(lo), hs[b]->FindBin(hi));
      double c = ctr[b];
      double dOm = 2*TMath::Pi()*(std::cos((c-dcm/2)*TMath::DegToRad()) - std::cos((c+dcm/2)*TMath::DegToRad()));
      double A = acc->GetBinContent(acc->FindBin(c));
      if (y <= 0 || A <= 0.05 || c < 18 || c > 148) continue;
      double d = y / A / dOm, e = std::sqrt(y) / A / dOm;
      X.push_back(c); Y.push_back(d); E.push_back(e);
      printf("  %7.1f %8.0f  %6.3f %12.4g %9.3g\n", c, y, A, d, e);
   }
   printf("  -> %zu measured points\n", X.size());

   // ---- confront every potential, one free scale each ---------------------------------------
   const int NP = 5;
   const char *pk[NP] = {"K", "V", "G", "P", "M"};
   const char *pn[NP] = {"KD03 (Koning-Delaroche)", "CH89 (Varner)", "Becchetti-Greenlees", "Perey", "Menet"};
   const char *pv[NP] = {"24<A<209: A=14 below", "4<A<209, 16<E<65: E below", "40<A: below",
                         "30<A<100: below", "40<A, 30<E<60: both below"};
   int col[NP] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
   TString pdir = here + "/../ptolemy/dat/";
   std::vector<TGraph *> G(NP);
   std::vector<double> SC(NP), RMS(NP), CHI(NP);

   printf("\n  ==== shape test over %.0f-%.0f deg, one free scale per potential ====\n", fitLo, fitHi);
   printf("  potential                 scale(=L)   rms ln(d/f)   chi2/N   validity for p+14C @ 11.6 MeV\n");
   for (int i = 0; i < NP; ++i) {
      G[i] = eo::rd(pdir + TString::Format("el_omp_%s.dat", pk[i]));
      // log-space scale: minimises rms of ln(data/fit), the scale-free shape metric
      double s = 0; int n = 0;
      for (size_t j = 0; j < X.size(); ++j) {
         if (X[j] < fitLo || X[j] > fitHi) continue;
         double f = G[i]->Eval(X[j]);
         if (f > 0 && Y[j] > 0) { s += std::log(Y[j] / f); ++n; }
      }
      double sc = std::exp(s / n);
      double r = 0, c2 = 0; int m = 0;
      for (size_t j = 0; j < X.size(); ++j) {
         if (X[j] < fitLo || X[j] > fitHi) continue;
         double f = sc * G[i]->Eval(X[j]);
         if (f <= 0 || Y[j] <= 0) continue;
         r += std::pow(std::log(Y[j] / f), 2);
         c2 += std::pow((Y[j] - f) / E[j], 2);
         ++m;
      }
      SC[i] = sc; RMS[i] = std::sqrt(r / m); CHI[i] = c2 / m;
      printf("  %-24s %10.2f %12.3f %9.1f   %s\n", pn[i], sc, RMS[i], CHI[i], pv[i]);
   }
   int best = 0;
   for (int i = 1; i < NP; ++i) if (RMS[i] < RMS[best]) best = i;
   printf("\n  best shape: %s (rms %.3f)\n", pn[best], RMS[best]);

   // ---- one panel per potential, as always ---------------------------------------------------
   auto *gd = new TGraphErrors(X.size(), &X[0], &Y[0], nullptr, &E[0]);
   gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->SetLineWidth(2);
   auto *c1 = new TCanvas("ceo", "", 1500, 900); c1->Divide(3, 2);
   for (int i = 0; i < NP; ++i) {
      c1->cd(i + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      auto *fr = gPad->DrawFrame(15, 0.3 * (*std::min_element(Y.begin(), Y.end())),
                                 152, 3.0 * (*std::max_element(Y.begin(), Y.end())));
      fr->SetTitle(Form("%s;#theta_{cm} [deg];d#sigma/d#Omega [arb.]", pn[i]));
      auto *q = new TGraph();
      for (int j = 0; j < G[i]->GetN(); ++j) q->SetPoint(j, G[i]->GetX()[j], SC[i] * G[i]->GetY()[j]);
      q->SetLineColor(col[i]); q->SetLineWidth(3); q->Draw("L same");
      gd->Draw("P same");
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.055);
      tx.DrawLatex(0.16, 0.26, Form("rms ln(d/f) = %.3f", RMS[i]));
      tx.DrawLatex(0.16, 0.19, Form("#chi^{2}/N = %.1f", CHI[i]));
      if (i == best) { tx.SetTextColor(kRed + 1); tx.DrawLatex(0.16, 0.12, "best shape"); }
   }
   // summary panel: all five, scaled, on the data
   c1->cd(6); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   auto *fr = gPad->DrawFrame(15, 0.3 * (*std::min_element(Y.begin(), Y.end())),
                              152, 3.0 * (*std::max_element(Y.begin(), Y.end())));
   fr->SetTitle("all five, each on its own best scale;#theta_{cm} [deg];d#sigma/d#Omega [arb.]");
   auto *lg = new TLegend(0.34, 0.63, 0.93, 0.89);
   lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
   for (int i = 0; i < NP; ++i) {
      auto *q = new TGraph();
      for (int j = 0; j < G[i]->GetN(); ++j) q->SetPoint(j, G[i]->GetX()[j], SC[i] * G[i]->GetY()[j]);
      q->SetLineColor(col[i]); q->SetLineWidth(2); q->Draw("L same");
      lg->AddEntry(q, Form("%s  (rms %.3f)", pn[i], RMS[i]), "l");
   }
   gd->Draw("P same"); lg->Draw();
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "04_elastic_data_vs_omp_" + tag + ".png");

   // save the measured distribution so this need not be recomputed
   TFile fo(here + "/plots/elastic_omp_" + tag + ".root", "RECREATE");
   gd->Write("elastic_measured");
   fo.Close();
   printf("  wrote %s04_elastic_data_vs_omp_%s.png\n", out.Data(), tag.Data());
}
