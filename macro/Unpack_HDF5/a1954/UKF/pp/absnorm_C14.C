/// @file absnorm_C14.C
/// @brief Put the excited-state angular distributions on an absolute scale, in mb/sr.
///
/// This is a normalisation and nothing else. It does not re-extract, re-fit or re-analyse the
/// elastic channel: it READS the existing g.s. result (elastic_sideband_C14.C, run yesterday,
/// left untouched) and uses it to measure the luminosity.
///
/// THE ONE PHYSICS INPUT. The excited states are confined to the first ~40 cm of the chamber
/// while the elastic fills all 95 cm. That is a trigger effect, not physics, and the simulation
/// does not contain the trigger -- so it is avoided rather than corrected, by comparing the two
/// over the same z window. The window 10-400 mm is the flat part of the z distribution, where the
/// trigger efficiency is constant and therefore cancels in the ratio.
///
/// Rather than re-running the elastic extraction inside that window, the existing g.s. yield is
/// SCALED by the fraction of elastic events that fall in it, measured per angular bin as a pure
/// count ratio with no fitting. The validated extraction is preserved exactly; only its
/// normalisation changes.
///
///     L(theta) = [ Y_gs(theta) * f(theta) / A_gs(theta) / dOmega(theta) ] / sigma_KD03(theta)
///
/// with f = N(z in window) / N(all z) for the elastic peak.
///
/// DEFORMATION LENGTHS. Every FRESCO inelastic curve was run with a reference deformation length
/// of 0.281 fm. The collective form factor is linear in delta, so the cross section goes as
/// delta^2 and a fitted scale k gives delta = 0.281 * sqrt(k). That scale is the only free
/// parameter left in the comparison, and it is a physical one.
///
///   root -b -q 'absnorm_C14.C()'

namespace
{
double dOmega(double loDeg, double hiDeg) // sr, integrated over phi
{
   return 2 * TMath::Pi() * (std::cos(loDeg * TMath::DegToRad()) - std::cos(hiDeg * TMath::DegToRad()));
}
TGraph *readFR(TString fn)
{
   auto *g = new TGraph();
   std::ifstream in(fn.Data());
   double a, x;
   int n = 0;
   while (in >> a >> x)
      g->SetPoint(n++, a, x);
   if (g->GetN() == 0) {
      printf("\033[1;31mcannot read %s\033[0m\n", fn.Data());
      return nullptr;
   }
   return g;
}
const double DELTA_REF = 0.281; // fm
} // namespace

void absnorm_C14(TString gsFile = "plots/elastic_sideband_gs_gf.root",
                 TString gsCache = "plots/proton_kin_300gfx_nc.root",
                 TString accGsDir = "/mnt/f/a1954_C14_acc_gf_z10_400/",
                 TString accExDir = "/mnt/f/a1954_C14_acc_gf_z10_400/", TString frDir = "",
                 Double_t zMin = 10.0, Double_t zMax = 400.0, Double_t gsLo = -0.5, Double_t gsHi = 0.7,
                 Double_t normLo = 75.0, Double_t normHi = 120.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frDir.IsNull())
      frDir = here + "/../fresco/outputs/";

   // ---- yesterday's g.s. yield, read, not recomputed ----
   TFile *fg = TFile::Open(here + "/" + gsFile);
   TH1D *Y = fg && !fg->IsZombie() ? (TH1D *)fg->Get("yield_sideband") : nullptr;
   if (!Y) {
      printf("\033[1;31mno yield_sideband in %s\033[0m\n", gsFile.Data());
      return;
   }
   Y = (TH1D *)Y->Clone("Ygs");
   Y->SetDirectory(nullptr);

   TFile *fc = TFile::Open(here + "/" + gsCache);
   TTree *t = fc && !fc->IsZombie() ? (TTree *)fc->Get("pk") : nullptr;
   TFile *fa = TFile::Open(accGsDir + "acceptance_merged_gs.root");
   TH1D *A = fa && !fa->IsZombie() ? (TH1D *)fa->Get("hAcc_gs_sum") : nullptr;
   TGraph *frEl = readFR(frDir + "p14C_el_161_dsdo.dat");
   if (!t || !A || !frEl) {
      printf("\033[1;31mmissing the cache, the g.s. acceptance or the KD03 elastic\033[0m\n");
      return;
   }

   printf("\n===== luminosity, from the EXISTING g.s. result scaled into z %.0f-%.0f mm =====\n", zMin, zMax);
   printf("  theta_cm |  Y_gs  |  f    | acc   |  KD03 [mb/sr] |  L [counts/mb]\n");
   std::vector<double> Lv;
   for (int b = 1; b <= Y->GetNbinsX(); ++b) {
      double lo = Y->GetBinLowEdge(b), w = Y->GetBinWidth(b), c = Y->GetBinCenter(b);
      double y = Y->GetBinContent(b);
      if (y <= 0)
         continue;
      TString base = TString::Format("ex>%g&&ex<%g&&thcm>=%g&&thcm<%g", gsLo, gsHi, lo, lo + w);
      double nAll = t->GetEntries(base);
      double nIn = t->GetEntries(base + TString::Format("&&vertexz>%g&&vertexz<%g", zMin, zMax));
      if (nAll < 10)
         continue;
      double f = nIn / nAll, a = A->GetBinContent(A->FindBin(c)), s = frEl->Eval(c);
      if (a <= 0.05 || s <= 0)
         continue;
      double L = y * f / a / dOmega(lo, lo + w) / s;
      printf("  %3.0f-%3.0f  | %6.0f | %.3f | %.3f | %13.4g | %14.2f\n", lo, lo + w, y, f, a, s, L);
      if (c >= normLo && c <= normHi)
         Lv.push_back(L);
   }
   if (Lv.empty()) {
      printf("\033[1;31mno elastic bins in the normalisation range\033[0m\n");
      return;
   }
   std::sort(Lv.begin(), Lv.end());
   const double L = Lv[Lv.size() / 2], Llo = Lv.front(), Lhi = Lv.back();
   printf("\n  L = %.2f counts/mb  (median over %.0f-%.0f deg, where the elastic acceptance is flat;\n"
          "      the %zu bins there span %.2f to %.2f = %+.0f%% / %+.0f%%)\n",
          L, normLo, normHi, Lv.size(), Llo, Lhi, 100 * (Llo / L - 1), 100 * (Lhi / L - 1));

   // ---- the excited states on that scale ----
   struct Lev {
      const char *name, *file, *fr, *accH, *accF;
      int idx, Lmult;
   };
   const std::vector<Lev> LV = {
      {"6.09 (1^{-})", "plots/exc_angdist_gfex.root", "p14C_inel_161_6094_L1_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 0, 1},
      {"6.70 (3^{-})", "plots/exc_angdist_gfex.root", "p14C_inel_161_6728_L3_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 1, 3},
      {"7.00 (2^{+})", "plots/exc_angdist_gfex.root", "p14C_inel_161_7012_L2_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 2, 2},
      {"8.53", "plots/exc8_angdist_hi.root", "p14C_inel_161_8317_L2_dsdo_ex2.dat", "hAcc_ex8_sum",
       "acceptance_merged_ex8.root", 0, 2},
      {"9.36", "plots/exc8_angdist_hi.root", "p14C_inel_161_9363_L2_dsdo_ex2.dat", "hAcc_ex8_sum",
       "acceptance_merged_ex8.root", 1, 2}};

   printf("\n===== absolute cross sections and deformation lengths =====\n");
   printf("  level        | L | sigma_int [mb] | delta [fm] | stat  |  norm      | chi2/n\n");
   TCanvas *c1 = new TCanvas("cab", "absolute", 1500, 900);
   c1->Divide(3, 2);
   int pad = 0;
   for (auto &lv : LV) {
      TFile *fy = TFile::Open(here + "/" + lv.file);
      TFile *fA = TFile::Open(accExDir + lv.accF);
      TH1D *hy = fy && !fy->IsZombie() ? (TH1D *)fy->Get(TString::Format("yield_%d", lv.idx)) : nullptr;
      TH1D *ac = fA && !fA->IsZombie() ? (TH1D *)fA->Get(lv.accH) : nullptr;
      TGraph *fr = readFR(frDir + lv.fr);
      if (!hy || !ac || !fr) {
         printf("  %-12s  missing input\n", lv.name);
         continue;
      }
      hy = (TH1D *)hy->Clone(TString::Format("y%d_%s", lv.idx, lv.accH));
      hy->SetDirectory(nullptr);

      auto *hx = (TH1D *)hy->Clone(TString::Format("x%d_%s", lv.idx, lv.accH));
      hx->Reset();
      hx->SetDirectory(nullptr);
      double sInt = 0;
      for (int b = 1; b <= hy->GetNbinsX(); ++b) {
         double lo = hy->GetBinLowEdge(b), w = hy->GetBinWidth(b), c = hy->GetBinCenter(b);
         double N = hy->GetBinContent(b), e = hy->GetBinError(b);
         if (N <= 0)
            continue;
         double a = ac->GetBinContent(ac->FindBin(c));
         if (a <= 0.05)
            continue;
         double dO = dOmega(lo, lo + w);
         hx->SetBinContent(b, N / a / dO / L);
         hx->SetBinError(b, e / a / dO / L);
         sInt += N / a / L;
      }

      double sn = 0, sd = 0, c2 = 0;
      int nd = 0;
      for (int b = 1; b <= hx->GetNbinsX(); ++b) {
         double y = hx->GetBinContent(b), e = hx->GetBinError(b);
         if (y <= 0 || e <= 0)
            continue;
         double f = fr->Eval(hx->GetBinCenter(b));
         if (f <= 0)
            continue;
         sn += y * f / (e * e);
         sd += f * f / (e * e);
      }
      double k = sd > 0 ? sn / sd : 0;
      for (int b = 1; b <= hx->GetNbinsX(); ++b) {
         double y = hx->GetBinContent(b), e = hx->GetBinError(b);
         if (y <= 0 || e <= 0)
            continue;
         double f = fr->Eval(hx->GetBinCenter(b)) * k;
         c2 += (y - f) * (y - f) / (e * e);
         ++nd;
      }
      double d = DELTA_REF * std::sqrt(std::max(k, 0.0));
      double eK = sd > 0 ? 1.0 / std::sqrt(sd) : 0;
      double ed = k > 0 ? 0.5 * d * eK / k : 0;
      printf("  %-12s | %d | %14.4g | %10.3f | %5.3f | %+.2f/%+.2f | %6.2f\n", lv.name, lv.Lmult, sInt, d, ed,
             d * (std::sqrt(L / Lhi) - 1), d * (std::sqrt(L / Llo) - 1), nd > 1 ? c2 / (nd - 1) : -1);

      if (pad < 6) {
         c1->cd(++pad);
         gPad->SetLogy();
         hx->SetTitle(TString::Format("%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", lv.name));
         hx->SetMarkerStyle(20);
         hx->SetMarkerSize(1.2);
         hx->SetLineWidth(2);
         hx->Draw("E1");
         auto *g2 = new TGraph(*fr);
         for (int n = 0; n < g2->GetN(); ++n)
            g2->SetPointY(n, g2->GetPointY(n) * k);
         g2->SetLineColor(kRed + 1);
         g2->SetLineWidth(3);
         g2->Draw("L same");
         auto *lg = new TLegend(0.42, 0.76, 0.89, 0.89);
         lg->AddEntry(hx, "data (absolute)", "lp");
         lg->AddEntry(g2, TString::Format("L=%d  #delta=%.2f fm", lv.Lmult, d), "l");
         lg->SetTextSize(0.042);
         lg->Draw();
      }
   }
   c1->SaveAs(here + "/plots/absnorm_C14.png");
   printf("\nwrote plots/absnorm_C14.png\n\n");
}
