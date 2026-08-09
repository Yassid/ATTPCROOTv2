/// @file absnorm_C14.C
/// @brief Put the excited-state angular distributions on an ABSOLUTE scale using the elastic.
///
/// This is a normalisation, not a re-analysis. The excited-state yields are taken exactly as
/// exc_angdist_C14.C and exc8_angdist_C14.C already extract them; all this does is measure the
/// scale that converts them to mb/sr, and then fit the deformation length with no free
/// normalisation left over.
///
/// THE ELASTIC YIELD COMES FROM THE SIDEBAND EXTRACTION, not from a fixed Ex window. The elastic
/// peak broadens with angle (sigma ~0.17 MeV forward to ~0.9 backward), so a fixed window slides
/// off it and the yield collapses; done that way the luminosity came out varying by a factor of
/// 9 across the angular range, which is not a luminosity but a measure of how much peak the
/// window was losing. elastic_sideband_C14.C tracks the locus and subtracts sidebands instead,
/// assuming no peak shape, and it is what reproduces the calculated elastic distribution.
///
/// HOW THE SCALE IS MEASURED. The luminosity (beam particles x target areal density) is not known
/// independently, but the elastic cross section is: the Koning-Delaroche global potential gives it
/// absolutely. So the elastic yield measures the luminosity, and the same luminosity is applied to
/// the inelastic yields.
///
/// WHY THE VERTEX-Z WINDOW MATTERS, and it is the whole reason this macro exists. The elastic
/// events fill the entire 95 cm of the chamber while the excited states are confined to the first
/// ~40 cm. That is a TRIGGER effect, not a physics one, and the simulation does not contain the
/// trigger -- so it cannot be corrected for, only avoided. It is avoided by counting both the
/// elastic and the inelastic in the SAME z window, chosen to be the flat part of the distribution
/// (10-400 mm) where the trigger efficiency is constant. In that window the flux, the target
/// thickness and the trigger efficiency are common to numerator and denominator and cancel in the
/// ratio. Integrating the elastic over 95 cm against an inelastic yield over 40 cm would instead
/// inflate every inelastic cross section by roughly the ratio of the two thicknesses.
///
/// The same window is where the vertex is measurable at all: in simulation the reconstructed z
/// tracks the truth to about 15 mm below 600 mm, but at a true z of 890 mm the mean reconstructed
/// z is 409 mm.
///
/// DEFORMATION LENGTHS. The FRESCO inelastic curves were computed with a reference deformation
/// length delta_ref = 0.281 fm. In the collective model the form factor is linear in delta, so the
/// cross section goes as delta^2 and a fitted scale k gives delta = delta_ref * sqrt(k). That is
/// the only free parameter left, and it is a physical one.
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
const double DELTA_REF = 0.281; // fm, the deformation length every inelastic curve was run with
} // namespace

void absnorm_C14(TString cache = "plots/proton_kin_300gfx_ex.root",
                 TString accDir = "/mnt/f/a1954_C14_acc_gf_z10_400/", TString frDir = "",
                 Double_t zMin = 10.0, Double_t zMax = 400.0, Double_t gsLo = -0.27, Double_t gsHi = 0.58,
                 Double_t normLo = 30.0, Double_t normHi = 55.0,
                 TString elFile = "plots/elastic_sideband_gs_absn.root")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frDir.IsNull())
      frDir = here + "/../fresco/outputs/";

   TFile *fd = TFile::Open(here + "/" + cache);
   TFile *fg = TFile::Open(accDir + "acceptance_merged_gs.root");
   if (!fd || fd->IsZombie() || !fg || fg->IsZombie()) {
      printf("\033[1;31mcannot open the cache or the g.s. acceptance\033[0m\n");
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *accGs = (TH1D *)fg->Get("hAcc_gs_sum");
   TGraph *frEl = readFR(frDir + "p14C_el_161_dsdo.dat");
   if (!t || !accGs || !frEl)
      return;

   // ================= 1. the luminosity, from the elastic in the z window =================
   TFile *fe = TFile::Open(here + "/" + elFile);
   TH1D *hEl = fe && !fe->IsZombie() ? (TH1D *)fe->Get("yield_sideband") : nullptr;
   if (!hEl) {
      printf("\033[1;31mno yield_sideband in %s -- run elastic_sideband_C14.C in this z window "
             "first\033[0m\n",
             elFile.Data());
      return;
   }
   hEl = (TH1D *)hEl->Clone("hEl");
   hEl->SetDirectory(nullptr);
   const int NE = hEl->GetNbinsX();

   printf("\n===== elastic, vertex z %.0f-%.0f mm =====\n", zMin, zMax);
   printf("  theta_cm |  counts | acc   | data/dOmega | KD03 [mb/sr] | L [counts/mb]\n");
   std::vector<double> Lv;
   for (int b = 1; b <= NE; ++b) {
      double lo = hEl->GetBinLowEdge(b), dEl = hEl->GetBinWidth(b), c = hEl->GetBinCenter(b);
      double N = hEl->GetBinContent(b), A = accGs->GetBinContent(accGs->FindBin(c));
      if (N < 5 || A <= 0.05)
         continue;
      double y = N / A / dOmega(lo, lo + dEl), s = frEl->Eval(c);
      if (s <= 0)
         continue;
      printf("  %3.0f-%3.0f  | %7.0f | %.3f | %11.4g | %12.4g | %12.2f\n", lo, lo + dEl, N, A, y, s, y / s);
      if (c >= normLo && c <= normHi)
         Lv.push_back(y / s);
   }
   if (Lv.empty()) {
      printf("\033[1;31mno elastic bins inside the normalisation range\033[0m\n");
      return;
   }
   std::sort(Lv.begin(), Lv.end());
   const double L = Lv[Lv.size() / 2];
   const double Llo = Lv.front(), Lhi = Lv.back();
   printf("\n  L = %.2f counts/mb   (median over %.0f-%.0f deg; the %zu bins there span %.2f to %.2f,\n"
          "      i.e. %+.0f%% / %+.0f%%, and that spread is the normalisation systematic)\n",
          L, normLo, normHi, Lv.size(), Llo, Lhi, 100 * (Llo / L - 1), 100 * (Lhi / L - 1));

   // ================= 2. the excited states, on that scale =================
   struct Lev {
      const char *name, *file, *frFile, *accName, *accFile;
      int idx, L;
   };
   // The multiplet uses the 6.094 acceptance (the only one measured for it); the two upper
   // structures use the 8.317 one where it exists. 7.27 is unnatural parity and has no curve.
   const std::vector<Lev> LV = {
      {"6.09 (1-)", "plots/exc_angdist_gfex.root", "p14C_inel_161_6094_L1_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 0, 1},
      {"6.70 (3-)", "plots/exc_angdist_gfex.root", "p14C_inel_161_6728_L3_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 1, 3},
      {"7.00 (2+)", "plots/exc_angdist_gfex.root", "p14C_inel_161_7012_L2_dsdo_ex2.dat", "hAcc_ex1_sum",
       "acceptance_merged_ex1.root", 2, 2},
      {"8.53", "plots/exc8_angdist_hi.root", "p14C_inel_161_8317_L2_dsdo_ex2.dat", "hAcc_ex8_sum",
       "acceptance_merged_ex8.root", 0, 2},
      {"9.36", "plots/exc8_angdist_hi.root", "p14C_inel_161_9363_L2_dsdo_ex2.dat", "hAcc_ex8_sum",
       "acceptance_merged_ex8.root", 1, 2}};

   printf("\n===== absolute cross sections and deformation lengths =====\n");
   printf("  level      | L | sigma_int [mb] | delta_L [fm] | stat  | norm syst | chi2/n\n");
   TCanvas *c1 = new TCanvas("cab", "absolute", 1500, 900);
   c1->Divide(3, 2);
   int pad = 0;

   for (auto &lv : LV) {
      TFile *fy = TFile::Open(here + "/" + lv.file);
      TFile *fa = TFile::Open(accDir + lv.accFile);
      if (!fy || fy->IsZombie() || !fa || fa->IsZombie()) {
         printf("  %-10s  missing %s or %s\n", lv.name, lv.file, lv.accFile);
         continue;
      }
      auto *hy = (TH1D *)fy->Get(TString::Format("yield_%d", lv.idx));
      auto *ac = (TH1D *)fa->Get(lv.accName);
      TGraph *fr = readFR(frDir + lv.frFile);
      if (!hy || !ac || !fr) {
         printf("  %-10s  missing yield/acceptance/curve\n", lv.name);
         continue;
      }
      hy = (TH1D *)hy->Clone(TString::Format("y_%s", lv.name));
      hy->SetDirectory(nullptr);

      // yield -> mb/sr
      auto *hx = (TH1D *)hy->Clone(TString::Format("x_%s", lv.name));
      hx->Reset();
      hx->SetDirectory(nullptr);
      double sInt = 0;
      for (int b = 1; b <= hy->GetNbinsX(); ++b) {
         double lo = hy->GetBinLowEdge(b), w = hy->GetBinWidth(b), c = hy->GetBinCenter(b);
         double N = hy->GetBinContent(b), e = hy->GetBinError(b);
         if (N <= 0)
            continue;
         double A = ac->GetBinContent(ac->FindBin(c));
         if (A <= 0.05)
            continue;
         double dO = dOmega(lo, lo + w);
         hx->SetBinContent(b, N / A / dO / L);
         hx->SetBinError(b, e / A / dO / L);
         sInt += N / A / L; // integrated over the bin's solid angle
      }

      // one free parameter: the scale of the calculation, i.e. (delta/delta_ref)^2
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
      double delta = DELTA_REF * std::sqrt(std::max(k, 0.0));
      double eK = sd > 0 ? 1.0 / std::sqrt(sd) : 0;             // error on the scale
      double eDelta = k > 0 ? 0.5 * delta * eK / k : 0;         // delta = ref*sqrt(k)
      // the luminosity spread propagates as sqrt, because sigma ~ 1/L and delta ~ sqrt(sigma)
      double dLo = delta * (std::sqrt(L / Lhi) - 1), dHi = delta * (std::sqrt(L / Llo) - 1);
      printf("  %-10s | %d | %14.4g | %12.4f | %5.4f | %+.3f/%+.3f | %6.2f\n", lv.name, lv.L, sInt, delta, eDelta,
             dLo, dHi, nd > 1 ? c2 / (nd - 1) : -1);

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
         auto *lg = new TLegend(0.45, 0.75, 0.89, 0.89);
         lg->AddEntry(hx, "data (absolute)", "lp");
         lg->AddEntry(g2, TString::Format("L=%d, #delta=%.3f fm", lv.L, delta), "l");
         lg->SetTextSize(0.04);
         lg->Draw();
      }
   }
   c1->SaveAs(here + "/plots/absnorm_C14.png");
   printf("\nwrote plots/absnorm_C14.png\n\n");
}
