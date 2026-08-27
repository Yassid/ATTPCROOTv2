/// @file mda_slices_C14.C
/// @brief Multipole decomposition analysis of the 14C(p,p') continuum, done the right way.
///
/// THE POINT OF AN MDA IS THAT NOTHING IS SEPARATED FIRST. The excitation-energy spectrum is cut
/// into slices; for each slice the RAW angular distribution is taken -- no peak fitting, no
/// assignment of counts to levels -- and decomposed as
///
///     dsigma/dOmega (theta ; Ex) = sum_L a_L(Ex) * sigma_L^DWBA(theta)
///
/// so the multipolarity comes out of the angular SHAPES and the 6.903 and 7.012 never have to be
/// resolved. Running a decomposition on per-level distributions that a Gaussian fit has already
/// separated presupposes exactly what the method exists to avoid; an earlier version of this
/// analysis did that and the result was meaningless.
///
/// The output is a_L against Ex: strength distributions per multipole, which is what shows whether
/// an L is present in a region without needing a level scheme.
///
/// THE BASIS. {L=2, L=3} only. Over 25-135 deg the normalised DWBA shapes correlate as L=0/L=2 =
/// 0.92 and L=2/L=4 = 0.96, so those three are one shape within the errors, while L=3 against L=2
/// is 0.30 and is separable. Condition numbers: {2,3} = 7.3 usable, {0,2,3,4} = 66 which returns
/// noise. So this asks whether octupole strength is present, and cannot search over all L.
///
/// SPANNING. Printed per slice: an MDA presupposes the basis can reach the data, and if the yield
/// falls far more steeply with angle than every basis shape then no positive combination can
/// represent it and the coefficients are an artefact of the basis being too flat.
///
///   root -b -q 'mda_slices_C14.C()'
namespace mds {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
struct Fit { double a2, a3, chi2; int ndf; };
Fit nnls(const std::vector<double> &y, const std::vector<double> &e,
         const std::vector<double> &f2, const std::vector<double> &f3)
{
   double s22 = 0, s33 = 0, s23 = 0, s2y = 0, s3y = 0;
   for (size_t i = 0; i < y.size(); ++i) {
      double w = 1.0 / (e[i] * e[i]);
      s22 += w * f2[i] * f2[i]; s33 += w * f3[i] * f3[i]; s23 += w * f2[i] * f3[i];
      s2y += w * f2[i] * y[i];  s3y += w * f3[i] * y[i];
   }
   double det = s22 * s33 - s23 * s23, a2 = 0, a3 = 0;
   if (std::fabs(det) > 1e-30) { a2 = (s33 * s2y - s23 * s3y) / det; a3 = (s22 * s3y - s23 * s2y) / det; }
   if (a2 < 0) { a2 = 0; a3 = s33 > 0 ? s3y / s33 : 0; }
   if (a3 < 0) { a3 = 0; a2 = s22 > 0 ? s2y / s22 : 0; }
   double c2 = 0;
   for (size_t i = 0; i < y.size(); ++i) c2 += std::pow((y[i] - a2 * f2[i] - a3 * f3[i]) / e[i], 2);
   return {a2, a3, c2, (int)y.size() - 2};
}
} // namespace mds

void mda_slices_C14(TString cache = "plots/proton_kin_cat5_s013.root",
                    TString accDir = "/mnt/f/a1954_C14_acc_catima_z10_490/",
                    Double_t exLo = 5.0, Double_t exHi = 10.5, Double_t dEx = 0.25,
                    Double_t cmLo = 25, Double_t cmHi = 135, Double_t dcm = 10,
                    Double_t vzLo = 10, Double_t vzHi = 490, Double_t chi2Cut = 5.0,
                    Double_t lumi = 72.5, Int_t nBoot = 400, UInt_t seed = 20260827)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TFile *fc = TFile::Open(here + "/" + cache);
   TNtuple *t = fc && !fc->IsZombie() ? (TNtuple *)fc->Get("pk") : nullptr;
   TFile *fa = TFile::Open(accDir + "acceptance_merged_ex1.root");
   TH1D *acc = fa && !fa->IsZombie() ? (TH1D *)fa->Get("hAcc_ex1_sum") : nullptr;
   if (!t || !acc) { printf("\033[1;31mmissing cache or acceptance\033[0m\n"); return; }

   const int NE = (int)std::lround((exHi - exLo) / dEx);
   const int NC = (int)std::lround((cmHi - cmLo) / dcm);
   // raw counts in (Ex, theta_cm). NOTHING is fitted or assigned here -- that is the point.
   auto *h = new TH2D("hmda", "raw yield;E_{x} [MeV];#theta_{cm} [deg]", NE, exLo, exHi, NC, cmLo, cmHi);
   float *v;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); v = t->GetArgs();
      if (v[5] > chi2Cut || v[2] < vzLo || v[2] > vzHi) continue;
      h->Fill(v[4], v[3]);
   }
   printf("\n  MDA on %0.f raw tracks -- no peak fitting, no level assignment\n", h->GetEntries());
   printf("  basis {L=2, L=3} at Ex = 7.012 MeV, Perey; slices of %.2f MeV over %.1f-%.1f\n",
          dEx, exLo, exHi);

   auto *g2 = mds::rd(pdir + "mdaP_7012_L2.dat");
   auto *g3 = mds::rd(pdir + "mdaP_7012_L3.dat");
   if (!g2->GetN() || !g3->GetN()) { printf("\033[1;31mno basis curves\033[0m\n"); return; }

   TRandom3 rng(seed);
   auto *hA2 = new TH1D("hA2", "L=2 strength;E_{x} [MeV];d#sigma/d#Omega integrated [mb]", NE, exLo, exHi);
   auto *hA3 = new TH1D("hA3", "L=3 strength", NE, exLo, exHi);
   printf("\n  %8s %10s %10s %10s %9s %8s  %s\n",
          "Ex [MeV]", "counts", "L=2 [mb]", "L=3 [mb]", "f(L=3)", "chi2/ndf", "spans?");
   for (int b = 1; b <= NE; ++b) {
      std::vector<double> y, e, f2, f3, th;
      double ntot = 0;
      for (int c = 1; c <= NC; ++c) {
         double cm = h->GetYaxis()->GetBinCenter(c);
         double n = h->GetBinContent(b, c);
         ntot += n;
         double A = acc->GetBinContent(acc->FindBin(cm));
         if (A <= 0.05) continue;
         double dOm = 2 * TMath::Pi() * (std::cos((cm - dcm / 2) * TMath::DegToRad())
                                       - std::cos((cm + dcm / 2) * TMath::DegToRad()));
         double xs = n / A / dOm / lumi;
         double er = std::sqrt(std::max(n, 1.0)) / A / dOm / lumi;   // sqrt(N), floored at 1
         th.push_back(cm); y.push_back(xs); e.push_back(er);
         f2.push_back(g2->Eval(cm)); f3.push_back(g3->Eval(cm));
      }
      if (th.size() < 4 || ntot < 20) continue;
      auto r = mds::nnls(y, e, f2, f3);
      double i2 = 0, i3 = 0;
      for (size_t k = 0; k < th.size(); ++k) {
         double w = 2 * TMath::Pi() * std::sin(th[k] * TMath::DegToRad()) * dcm * TMath::DegToRad();
         i2 += r.a2 * f2[k] * w; i3 += r.a3 * f3[k] * w;
      }
      std::vector<double> bf;
      for (int q = 0; q < nBoot; ++q) {
         std::vector<double> yy(y.size());
         for (size_t k = 0; k < y.size(); ++k) yy[k] = rng.Gaus(y[k], e[k]);
         auto rr = mds::nnls(yy, e, f2, f3);
         double j2 = 0, j3 = 0;
         for (size_t k = 0; k < th.size(); ++k) { j2 += rr.a2 * f2[k]; j3 += rr.a3 * f3[k]; }
         if (j2 + j3 > 0) bf.push_back(j3 / (j2 + j3));
      }
      double m = 0, s = 0;
      for (double u : bf) m += u; m /= std::max<size_t>(1, bf.size());
      for (double u : bf) s += (u - m) * (u - m);
      s = std::sqrt(s / std::max<size_t>(1, bf.size()));
      double dFall = y.front() / std::max(1e-12, y.back());
      double f2F = f2.front() / f2.back(), f3F = f3.front() / f3.back();
      bool spans = dFall <= 1.5 * std::max(f2F, f3F) && dFall >= 0.67 * std::min(f2F, f3F);
      hA2->SetBinContent(b, i2); hA3->SetBinContent(b, i3);
      printf("  %8.2f %10.0f %10.3f %10.3f %5.2f+-%.2f %8.2f  %s\n",
             h->GetXaxis()->GetBinCenter(b), ntot, i2, i3, i3 / std::max(1e-12, i2 + i3), s,
             r.ndf > 0 ? r.chi2 / r.ndf : 0, spans ? "yes" : "NO");
   }

   auto *c1 = new TCanvas("cmds", "", 1300, 560); c1->Divide(2, 1);
   c1->cd(1); gPad->SetRightMargin(0.13); gPad->SetLogz(); h->Draw("colz");
   c1->cd(2); gPad->SetGridx(); gPad->SetGridy();
   hA2->SetLineColor(kAzure + 2); hA2->SetLineWidth(3); hA2->SetFillColorAlpha(kAzure + 2, 0.25);
   hA3->SetLineColor(kGreen + 3); hA3->SetLineWidth(3);
   hA2->SetTitle("multipole strength vs excitation energy;E_{x} [MeV];#sigma integrated [mb]");
   hA2->Draw("hist"); hA3->Draw("hist same");
   auto *lg = new TLegend(0.55, 0.72, 0.90, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hA2, "L = 2", "l"); lg->AddEntry(hA3, "L = 3", "l"); lg->Draw();
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "18_mda_slices.png");
   printf("\n  wrote %s18_mda_slices.png\n\n", out.Data());
}
