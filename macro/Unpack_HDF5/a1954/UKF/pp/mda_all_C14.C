/// @file mda_all_C14.C
/// @brief MDA of the 6.3-7.6 MeV region with the FULL multipole basis, one panel per slice.
///
/// The previous pass restricted the basis to {L=2,3} on conditioning grounds and never asked the
/// prior question: can ANY non-negative combination of multipoles reproduce the measured angular
/// shape? That question is well posed even when the individual coefficients are not, because the
/// FIT QUALITY of the best combination does not care about degeneracy -- only the decomposition
/// into components does. So here the basis is L = 0,1,2,3,4,5 and the headline number per slice is
/// chi2/ndf of the best non-negative combination.
///
/// Every subset of the basis is also fitted exhaustively (63 of them), which says how many
/// multipoles the shape actually demands rather than how many it tolerates.
///
///   root -b -q 'mda_all_C14.C()'
namespace mall {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
/// Non-negative least squares by coordinate descent. Converges for any number of components;
/// unlike the closed-form 2x2 it does not need a case analysis on which coefficient goes negative.
std::vector<double> nnls(const std::vector<std::vector<double>> &F, const std::vector<double> &y,
                         const std::vector<double> &e)
{
   const size_t M = F.size(), N = y.size();
   std::vector<double> a(M, 0.0), sjj(M, 0.0);
   std::vector<std::vector<double>> sjk(M, std::vector<double>(M, 0.0));
   std::vector<double> sjy(M, 0.0);
   for (size_t j = 0; j < M; ++j) {
      for (size_t i = 0; i < N; ++i) {
         double w = 1.0 / (e[i] * e[i]);
         sjy[j] += w * F[j][i] * y[i];
         for (size_t k = 0; k < M; ++k) sjk[j][k] += w * F[j][i] * F[k][i];
      }
      sjj[j] = sjk[j][j];
   }
   for (int it = 0; it < 2000; ++it) {
      double mv = 0;
      for (size_t j = 0; j < M; ++j) {
         if (sjj[j] <= 0) continue;
         double num = sjy[j];
         for (size_t k = 0; k < M; ++k) if (k != j) num -= a[k] * sjk[j][k];
         double nv = std::max(0.0, num / sjj[j]);
         mv = std::max(mv, std::fabs(nv - a[j]));
         a[j] = nv;
      }
      if (mv < 1e-12) break;
   }
   return a;
}
double chi2of(const std::vector<std::vector<double>> &F, const std::vector<double> &a,
              const std::vector<double> &y, const std::vector<double> &e)
{
   double c = 0;
   for (size_t i = 0; i < y.size(); ++i) {
      double m = 0;
      for (size_t j = 0; j < F.size(); ++j) m += a[j] * F[j][i];
      c += std::pow((y[i] - m) / e[i], 2);
   }
   return c;
}
} // namespace mall

void mda_all_C14(TString cache = "plots/proton_kin_cat5_s013.root",
                 TString accDir = "/mnt/f/a1954_C14_acc_catima_z10_490/",
                 Double_t exLo = 6.3, Double_t exHi = 7.7, Double_t dEx = 0.20,
                 Double_t cmLo = 25, Double_t cmHi = 135, Double_t dcm = 10,
                 Double_t vzLo = 10, Double_t vzHi = 490, Double_t chi2Cut = 5.0,
                 Double_t lumi = 72.5, TString use = "0123",
                 // which set of basis curves: mdaAll = Perey, mdaV = CH89. The multipole answer MUST
                 // be checked against more than one potential -- see the retracted swap finding.
                 TString basis = "mdaAll",
                 Int_t nBoot = 500, UInt_t seed = 20260828)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TFile *fc = TFile::Open(here + "/" + cache);
   TNtuple *t = fc && !fc->IsZombie() ? (TNtuple *)fc->Get("pk") : nullptr;
   TFile *fa = TFile::Open(accDir + "acceptance_merged_ex1.root");
   TH1D *acc = fa && !fa->IsZombie() ? (TH1D *)fa->Get("hAcc_ex1_sum") : nullptr;
   if (!t || !acc) { printf("\033[1;31mmissing cache or acceptance\033[0m\n"); return; }

   // Only the multipoles named in `use` enter the basis. Over 25-135 deg the normalised shapes
   // correlate as L0/L4 = 0.95, L2/L4 = 0.94 and L3/L5 = 0.96, so L4 and L5 carry no shape that
   // L0/L2 and L3 do not already carry; including them only splits one shape's strength across
   // two indistinguishable labels. The default {0,1,2,3} is the largest basis whose members are
   // actually separable by these data.
   std::vector<TGraph *> gL; std::vector<int> Lid;
   const int colAll[6] = {kGray + 2, kOrange + 7, kAzure + 2, kGreen + 3, kMagenta + 2, kCyan + 2};
   std::vector<int> col;
   for (int L = 0; L < 6; ++L) {
      if (!use.Contains(Form("%d", L))) continue;
      auto *g = mall::rd(Form("%s%s_L%d.dat", pdir.Data(), basis.Data(), L));
      if (!g->GetN()) { printf("\033[1;31mmissing basis L=%d\033[0m\n", L); return; }
      gL.push_back(g); Lid.push_back(L); col.push_back(colAll[L]);
   }
   const int NL = (int)gL.size();

   const int NE = (int)std::lround((exHi - exLo) / dEx);
   const int NC = (int)std::lround((cmHi - cmLo) / dcm);
   auto *h = new TH2D("hall", "raw yield;E_{x} [MeV];#theta_{cm} [deg]", NE, exLo, exHi, NC, cmLo, cmHi);
   float *v;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); v = t->GetArgs();
      if (v[5] > chi2Cut || v[2] < vzLo || v[2] > vzHi) continue;
      h->Fill(v[4], v[3]);
   }
   printf("\n  MDA with basis L = {%s}, Perey, slices of %.0f keV over %.1f-%.1f MeV\n",
          use.Data(), dEx * 1000, exLo, exHi);
   printf("  %0.f tracks; no peak fitting, no level assignment\n\n", h->GetEntries());

   auto *cv = new TCanvas("call", "", 1500, 900);
   cv->Divide(3, 3);
   printf("  strength per multipole, integrated over the fitted angular range [mb].\n");
   printf("  errors are 16/84 bootstrap percentiles (%d resamples); a zero with a one-sided error\n", nBoot);
   printf("  means the NNLS boundary -- consistent with no strength, not a measurement of none.\n\n");
   printf("  %8s %7s |", "Ex", "counts");
   for (int L = 0; L < NL; ++L) printf(" %16s |", Form("L = %d", Lid[L]));
   printf(" %8s\n", "chi2/ndf");
   printf("  %s\n", TString('-', 34 + 19 * NL).Data());

   TRandom3 rng(seed);
   std::vector<TH1D *> hS(NL);
   for (int L = 0; L < NL; ++L) {
      hS[L] = new TH1D(Form("hS%d", Lid[L]), Form("L = %d;E_{x} [MeV];#sigma integrated [mb]", Lid[L]),
                       NE, exLo, exHi);
      hS[L]->SetLineColor(col[L]); hS[L]->SetLineWidth(3); hS[L]->SetMarkerColor(col[L]);
   }
   std::vector<double> tot(NL, 0), totLo(NL, 0), totHi(NL, 0);

   for (int b = 1; b <= NE; ++b) {
      std::vector<double> y, e, th;
      double ntot = 0;
      for (int c = 1; c <= NC; ++c) {
         double cm = h->GetYaxis()->GetBinCenter(c);
         double n = h->GetBinContent(b, c);
         ntot += n;
         double A = acc->GetBinContent(acc->FindBin(cm));
         if (A <= 0.05) continue;
         double dOm = 2 * TMath::Pi() * (std::cos((cm - dcm / 2) * TMath::DegToRad())
                                       - std::cos((cm + dcm / 2) * TMath::DegToRad()));
         th.push_back(cm);
         y.push_back(n / A / dOm / lumi);
         e.push_back(std::sqrt(std::max(n, 1.0)) / A / dOm / lumi);
      }
      if (th.size() < 4 || ntot < 20) continue;
      std::vector<std::vector<double>> F(NL, std::vector<double>(th.size()));
      for (int L = 0; L < NL; ++L)
         for (size_t k = 0; k < th.size(); ++k) F[L][k] = gL[L]->Eval(th[k]);

      // solid angle of each theta bin, for integrating a component to a cross section
      std::vector<double> w(th.size());
      for (size_t k = 0; k < th.size(); ++k)
         w[k] = 2 * TMath::Pi() * (std::cos((th[k] - dcm / 2) * TMath::DegToRad())
                                 - std::cos((th[k] + dcm / 2) * TMath::DegToRad()));
      auto integ = [&](const std::vector<double> &a, int L) {
         double s = 0; for (size_t k = 0; k < th.size(); ++k) s += a[L] * F[L][k] * w[k]; return s; };

      auto aFit = mall::nnls(F, y, e);
      int nUsed = 0; for (int L = 0; L < NL; ++L) if (aFit[L] > 0) ++nUsed;
      double c2 = mall::chi2of(F, aFit, y, e);
      int ndf = (int)y.size() - std::max(1, nUsed);

      // bootstrap: resample the data within errors, refit, keep the integrated strength per L
      std::vector<std::vector<double>> bs(NL);
      for (int q = 0; q < nBoot; ++q) {
         std::vector<double> yy(y.size());
         for (size_t k = 0; k < y.size(); ++k) yy[k] = rng.Gaus(y[k], e[k]);
         auto ab = mall::nnls(F, yy, e);
         for (int L = 0; L < NL; ++L) bs[L].push_back(integ(ab, L));
      }
      printf("  %8.2f %7.0f |", h->GetXaxis()->GetBinCenter(b), ntot);
      for (int L = 0; L < NL; ++L) {
         std::sort(bs[L].begin(), bs[L].end());
         double lo = bs[L][(size_t)(0.16 * bs[L].size())], hi = bs[L][(size_t)(0.84 * bs[L].size())];
         double v0 = integ(aFit, L);
         hS[L]->SetBinContent(b, v0);
         hS[L]->SetBinError(b, 0.5 * (hi - lo));
         tot[L] += v0; totLo[L] += lo; totHi[L] += hi;
         printf(" %6.2f +%5.2f-%5.2f |", v0, std::max(0.0, hi - v0), std::max(0.0, v0 - lo));
      }
      printf(" %8.2f\n", ndf > 0 ? c2 / ndf : 0);

      cv->cd(b); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      auto *gd = new TGraphErrors((int)th.size());
      for (size_t k = 0; k < th.size(); ++k) { gd->SetPoint(k, th[k], y[k]); gd->SetPointError(k, 0, e[k]); }
      gd->SetMarkerStyle(20); gd->SetMarkerSize(0.9);
      gd->SetTitle(Form("E_{x} = %.2f MeV  (%.0f counts);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]",
                        h->GetXaxis()->GetBinCenter(b), ntot));
      gd->Draw("AP");
      auto *gt = new TGraph((int)th.size());
      for (size_t k = 0; k < th.size(); ++k) {
         double m = 0; for (int L = 0; L < NL; ++L) m += aFit[L] * F[L][k];
         gt->SetPoint(k, th[k], m);
      }
      gt->SetLineColor(kRed + 1); gt->SetLineWidth(3); gt->Draw("L same");
      auto *lg = new TLegend(0.58, 0.62, 0.93, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->SetTextSize(0.035);
      lg->AddEntry(gt, Form("sum  #chi^{2}/ndf %.1f", ndf > 0 ? c2 / ndf : 0), "l");
      for (int L = 0; L < NL; ++L) {
         if (aFit[L] <= 0) continue;
         auto *gc = new TGraph((int)th.size());
         for (size_t k = 0; k < th.size(); ++k) gc->SetPoint(k, th[k], aFit[L] * F[L][k]);
         gc->SetLineColor(col[L]); gc->SetLineStyle(2); gc->SetLineWidth(2); gc->Draw("L same");
         lg->AddEntry(gc, Form("L = %d", Lid[L]), "l");
      }
      lg->Draw();
   }
   printf("  %s\n", TString('-', 34 + 19 * NL).Data());
   double sT = 0; for (int L = 0; L < NL; ++L) sT += tot[L];
   printf("  %8s %7s |", "SUM", "");
   for (int L = 0; L < NL; ++L) printf(" %6.2f +%5.2f-%5.2f |", tot[L],
          std::max(0.0, totHi[L] - tot[L]), std::max(0.0, tot[L] - totLo[L]));
   printf("\n  %8s %7s |", "fraction", "");
   for (int L = 0; L < NL; ++L) printf(" %15.3f  |", sT > 0 ? tot[L] / sT : 0);
   printf("\n\n");

   // strength distributions, one panel per multipole
   auto *cs = new TCanvas("cstr", "", 1300, 900); cs->Divide(2, 2);
   for (int L = 0; L < NL; ++L) {
      cs->cd(L + 1); gPad->SetGridx(); gPad->SetGridy();
      hS[L]->SetMarkerStyle(20);
      hS[L]->SetTitle(Form("L = %d strength;E_{x} [MeV];#sigma integrated [mb]", Lid[L]));
      hS[L]->Draw("E1");
   }
   cs->SaveAs("/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/22_mda_strength_perL.png");

   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   cv->SaveAs(out + "19_mda_allL_slices.png");
   printf("\n  wrote %s19_mda_allL_slices.png\n\n", out.Data());
}
