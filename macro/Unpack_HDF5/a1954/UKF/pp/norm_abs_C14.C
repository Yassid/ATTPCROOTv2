/// @file norm_abs_C14.C
/// @brief Fix the ABSOLUTE normalisation of the a1954 14C(p,p') data from the elastic channel.
///
/// Everything downstream of this macro is quoted in mb/sr rather than in arbitrary units, and
/// every excited-state deformation length depends on the single constant it produces.
///
/// THE IDEA. The luminosity (beam particles x target areal density) is not known independently
/// here -- there is no separate beam counter measurement to rely on. But the elastic cross
/// section IS known: the Koning-Delaroche global optical potential predicts it absolutely, and
/// the same calculation already reproduces the measured elastic shape including the diffraction
/// minimum. So the elastic yield is used to measure the luminosity, and that luminosity is then
/// applied unchanged to the inelastic yields. What comes out is a ratio to the elastic cross
/// section, with all its systematics, but it is an absolute number.
///
/// THE VERTEX-Z SLAB IS THE POINT OF THIS MACRO, and it is easy to get silently wrong. The
/// elastic events fill the whole 95 cm of the chamber; the excited states populate only the first
/// ~40 cm and are gone by 60 cm. Taking the ratio of an elastic yield integrated over 95 cm to an
/// inelastic yield integrated over 40 cm compares two different target thicknesses and would
/// inflate every inelastic cross section by roughly the ratio of the two. Both are therefore
/// restricted to the SAME z window, so that the luminosity cancels exactly rather than
/// approximately, and the acceptance is recomputed in that same window
/// (diagnostics/acc_zslab.sh) so the correction describes the same selection.
///
/// The window 10-400 mm is also where the vertex is measurable: reconstructed z tracks the truth
/// to about 15 mm below 600 mm, but at a true z of 890 mm the mean reconstructed z is 409 mm.
///
/// UNITS. Counts are histogrammed in equal theta_cm bins, so the solid angle of a bin is
/// dOmega = 2*pi*(cos(theta_lo) - cos(theta_hi)) steradian, NOT a constant. Dividing by it is
/// what turns a yield into a cross section; forgetting it tilts the whole distribution by
/// sin(theta).
///
///   root -b -q 'norm_abs_C14.C()'

namespace
{
// solid angle of a theta bin, integrated over phi
double dOmega(double loDeg, double hiDeg)
{
   return 2 * TMath::Pi() * (std::cos(loDeg * TMath::DegToRad()) - std::cos(hiDeg * TMath::DegToRad()));
}
} // namespace

void norm_abs_C14(TString cache = "plots/proton_kin_300gfx_ex.root",
                  TString accDir = "/mnt/f/a1954_C14_acc_gf_z10_400/", TString frDir = "",
                  Double_t zMin = 10.0, Double_t zMax = 400.0, Double_t gsLo = -0.27, Double_t gsHi = 0.58,
                  Double_t cmMin = 20.0, Double_t cmMax = 150.0, Double_t dcm = 5.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frDir.IsNull())
      frDir = here + "/../fresco/outputs/";

   TFile *fd = TFile::Open(here + "/" + cache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   if (!fd || fd->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mcannot open the cache or the g.s. acceptance in %s\033[0m\n", accDir.Data());
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   if (!t || !acc) {
      printf("\033[1;31mmissing pk or hAcc_gs_sum\033[0m\n");
      return;
   }

   // ---- the elastic yield in the slab ----
   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   auto *hY = new TH1D("hY", "elastic yield;#theta_{cm} [deg];counts", NB, cmMin, cmMax);
   hY->Sumw2();
   t->Draw("thcm>>hY", TString::Format("ex>%g&&ex<%g&&vertexz>%g&&vertexz<%g", gsLo, gsHi, zMin, zMax), "goff");
   hY->SetDirectory(nullptr);

   // ---- FRESCO elastic, absolute mb/sr ----
   auto *fr = new TGraph();
   {
      std::ifstream in((frDir + "p14C_el_161_dsdo.dat").Data());
      double a, x;
      int n = 0;
      while (in >> a >> x)
         fr->SetPoint(n++, a, x);
      if (fr->GetN() == 0) {
         printf("\033[1;31mcannot read the FRESCO elastic\033[0m\n");
         return;
      }
   }

   // ---- yield -> arbitrary-unit cross section, then the luminosity ----
   auto *hD = new TH1D("hD", "d#sigma/d#Omega;#theta_{cm} [deg];arb", NB, cmMin, cmMax);
   hD->Sumw2();
   printf("\n===== elastic in the vertex-z slab %.0f-%.0f mm =====\n", zMin, zMax);
   printf("  theta_cm |  counts | acc   | dOmega [sr] | data (arb) | FRESCO [mb/sr] | ratio\n");
   for (int b = 1; b <= NB; ++b) {
      double lo = hY->GetBinLowEdge(b), hi = lo + dcm, c = hY->GetBinCenter(b);
      double N = hY->GetBinContent(b), eN = hY->GetBinError(b);
      double A = acc->GetBinContent(acc->FindBin(c));
      double dO = dOmega(lo, hi);
      if (N <= 0 || A <= 0.05 || dO <= 0)
         continue;
      double y = N / A / dO, ey = eN / A / dO;
      hD->SetBinContent(b, y);
      hD->SetBinError(b, ey);
      double th = fr->Eval(c);
      printf("  %3.0f-%3.0f  | %7.0f | %.3f | %11.4f | %10.3g | %14.4g | %8.4g\n", lo, hi, N, A, dO, y, th,
             th > 0 ? y / th : 0);
   }

   // ---- the luminosity, over several angular ranges ----
   // One number from one range would hide the fact that it is range-dependent, which is the
   // dominant systematic on every absolute cross section derived from it.
   struct Rng {
      double lo, hi;
      const char *why;
   };
   const std::vector<Rng> RS = {{25, 50, "forward: largest cross section, before the minimum"},
                                {30, 50, "forward, avoiding the lowest-acceptance bins"},
                                {70, 110, "the plateau, away from the diffraction minimum"},
                                {70, 130, "plateau, extended backward"},
                                {25, 130, "everything usable"}};
   printf("\n===== luminosity from each angular range =====\n");
   printf("  range    |   L [counts/mb] | rms of ln(data/L*FRESCO) | n | note\n");
   std::vector<double> Ls;
   for (auto &r : RS) {
      double sn = 0, sd = 0;
      int n = 0;
      for (int b = 1; b <= NB; ++b) {
         double c = hD->GetBinCenter(b), y = hD->GetBinContent(b), e = hD->GetBinError(b);
         if (c < r.lo || c > r.hi || y <= 0 || e <= 0)
            continue;
         double f = fr->Eval(c);
         if (f <= 0)
            continue;
         sn += y * f / (e * e);
         sd += f * f / (e * e);
         ++n;
      }
      if (sd <= 0)
         continue;
      double L = sn / sd, s2 = 0;
      int m = 0;
      for (int b = 1; b <= NB; ++b) {
         double c = hD->GetBinCenter(b), y = hD->GetBinContent(b);
         if (c < r.lo || c > r.hi || y <= 0)
            continue;
         double f = fr->Eval(c) * L;
         if (f <= 0)
            continue;
         s2 += std::pow(std::log(y / f), 2);
         ++m;
      }
      Ls.push_back(L);
      printf("  %3.0f-%3.0f  | %15.4g | %24.3f | %d | %s\n", r.lo, r.hi, L, m ? std::sqrt(s2 / m) : 0, n, r.why);
   }

   if (Ls.size() > 1) {
      double mn = *std::min_element(Ls.begin(), Ls.end()), mx = *std::max_element(Ls.begin(), Ls.end());
      double md = 0;
      { std::vector<double> v = Ls; std::sort(v.begin(), v.end()); md = v[v.size() / 2]; }
      printf("\n  median L = %.4g counts/mb; the ranges span %.4g to %.4g, i.e. %+.0f%% / %+.0f%%.\n", md, mn, mx,
             100 * (mn / md - 1), 100 * (mx / md - 1));
      printf("  \033[1mThat spread is the normalisation systematic\033[0m and it propagates directly into every\n"
             "  deformation length as its square root.\n");
   }

   // ---- store it ----
   TFile fo(here + "/plots/norm_abs_C14.root", "RECREATE");
   hY->Write("yield_gs");
   hD->Write("dsdo_gs_arb");
   if (!Ls.empty()) {
      std::vector<double> v = Ls;
      std::sort(v.begin(), v.end());
      TVectorD lum(3);
      lum[0] = v[v.size() / 2];
      lum[1] = v.front();
      lum[2] = v.back();
      lum.Write("luminosity_median_min_max");
   }
   TNamed slab(TString("zslab"), TString::Format("%.0f,%.0f", zMin, zMax));
   slab.Write();
   fo.Close();

   // ---- picture ----
   TCanvas *c1 = new TCanvas("cn", "absolute normalisation", 900, 620);
   gPad->SetLogy();
   double L = Ls.empty() ? 1 : [&] { std::vector<double> v = Ls; std::sort(v.begin(), v.end()); return v[v.size() / 2]; }();
   auto *hAbs = (TH1D *)hD->Clone("hAbs");
   hAbs->Scale(1.0 / L);
   hAbs->SetTitle(TString::Format("elastic, absolutely normalised (z %.0f-%.0f mm);#theta_{cm} [deg];"
                                  "d#sigma/d#Omega [mb/sr]",
                                  zMin, zMax));
   hAbs->SetMarkerStyle(20);
   hAbs->SetMarkerSize(1.1);
   hAbs->SetLineWidth(2);
   hAbs->Draw("E1");
   fr->SetLineColor(kRed + 1);
   fr->SetLineWidth(3);
   fr->Draw("L same");
   auto *lg = new TLegend(0.6, 0.75, 0.89, 0.88);
   lg->AddEntry(hAbs, "data", "lp");
   lg->AddEntry(fr, "KD03 elastic", "l");
   lg->Draw();
   c1->SaveAs(here + "/plots/norm_abs_C14.png");
   printf("\nwrote plots/norm_abs_C14.root and .png\n\n");
}
