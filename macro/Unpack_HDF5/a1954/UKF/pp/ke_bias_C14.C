/// @file ke_bias_C14.C
/// @brief Measure the reconstructed proton KE bias against elastic two-body kinematics.
///
/// For 14C(p,p) elastic the recoil proton energy is fixed by its lab angle alone: given the beam
/// energy, KE_elastic(theta_lab) is a single number. So every elastic event is a calibration
/// point, and the E_x drift seen in ex_slices_C14.C can be converted into the quantity that
/// actually causes it -- the error on the reconstructed kinetic energy as a function of energy.
///
/// This is more useful than E_x vs theta because it is stated in the fitter's own variable:
///   * a constant relative offset points at a scale error (B field, drift velocity, KE calibration)
///   * a bias growing with KE points at the energy-loss / range correction, which is applied over
///     the track length and so grows with how far the proton travels
///   * a bias flat in KE but varying with theta_lab points at the angle reconstruction
///
/// Method: per theta_lab bin, take events within a band of the running E_x locus (so the sample
/// is elastic, not the inelastic group), and compare their KE to KE_elastic(theta_lab), which is
/// obtained by bisecting the SAME two-body expression the analysis uses until E_x = 0.
///
///   root -b -q 'ke_bias_C14.C()'

#include <tuple>

static double kb_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static std::tuple<double, double> kb_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * kb_omega2(s, m1 * m1, m2 * m2) * kb_omega2(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg <= 0)
      return {NAN, NAN};
   double m4_ex = std::sqrt(arg);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double ct = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
               (kb_omega2(s, m1 * m1, m2 * m2) * kb_omega2(s, m3 * m3, m4_ex * m4_ex));
   ct = std::max(-1.0, std::min(1.0, ct));
   return {Ex, (TMath::Pi() - std::acos(ct)) * TMath::RadToDeg()};
}

/// KE of the elastic recoil proton at this lab angle: bisect until Ex = 0
static double keElastic(double m1, double m2, double Ebeam, double thlab)
{
   double lo = 0.01, hi = 0.999 * Ebeam;
   auto ex = [&](double k) { return std::get<0>(kb_kine(m1, m2, m2, m1, Ebeam, thlab, k)); };
   double flo = ex(lo), fhi = ex(hi);
   if (!std::isfinite(flo) || !std::isfinite(fhi))
      return NAN;
   if (flo * fhi > 0)
      return NAN;
   for (int i = 0; i < 200; ++i) {
      double mid = 0.5 * (lo + hi), fm = ex(mid);
      if (!std::isfinite(fm))
         return NAN;
      if (flo * fm <= 0) {
         hi = mid;
         fhi = fm;
      } else {
         lo = mid;
         flo = fm;
      }
   }
   return 0.5 * (lo + hi);
}

void ke_bias_C14(Double_t Ebeam = 161.0, Double_t thLo = 20.0, Double_t thHi = 88.0, Double_t dth = 4.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;

   const int NS = 2;
   const char *file[NS] = {"plots/proton_kin_300_ukf_nc.root", "plots/proton_kin_300gfx_nc.root"};
   const char *lbl[NS] = {"UKF", "GENFIT"};
   const int col[NS] = {kAzure + 2, kRed + 1};

   TFile *f[NS];
   TTree *t[NS];
   for (int i = 0; i < NS; ++i) {
      f[i] = TFile::Open(here + "/" + file[i]);
      if (!f[i] || f[i]->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", file[i]);
         return;
      }
      t[i] = (TTree *)f[i]->Get("pk");
   }

   const int NB = (int)((thHi - thLo) / dth);
   auto *gKE = new TGraphErrors[NS], *gRel = new TGraphErrors[NS], *gVsKE = new TGraphErrors[NS];
   int np[NS] = {0, 0};

   printf("\n===== KE bias vs theta_lab, elastic band, Ebeam = %.1f MeV =====\n", Ebeam);
   printf("  theta_lab | KE_elastic |        UKF: <KE>   bias   rel    N |     GENFIT: <KE>   bias   rel    N\n");
   for (int b = 0; b < NB; ++b) {
      double lo = thLo + b * dth, hi = lo + dth, ctr = 0.5 * (lo + hi);
      double keEl = keElastic(m_C14, m_p, Ebeam, ctr * TMath::DegToRad());
      if (!std::isfinite(keEl))
         continue;
      printf("  %4.0f-%-4.0f | %10.3f |", lo, hi, keEl);
      for (int i = 0; i < NS; ++i) {
         // elastic band: the KE that would give Ex=0, +-12 % (wide enough to hold the drifted
         // peak, narrow enough to exclude the 6-7 MeV group which sits ~25 % lower in KE)
         TString cut = TString::Format("theta>=%g&&theta<%g&&ke>%g&&ke<%g", lo, hi, 0.88 * keEl, 1.12 * keEl);
         auto *h = new TH1D(TString::Format("hk%d_%d", i, b), "", 120, 0.85 * keEl, 1.15 * keEl);
         t[i]->Draw(TString::Format("ke>>hk%d_%d", i, b), cut, "goff");
         h->SetDirectory(nullptr);
         double N = h->Integral();
         if (N < 40) {
            printf(" %38s |", "too few");
            delete h;
            continue;
         }
         // gaussian on the elastic line
         int bm = h->GetMaximumBin();
         double c0 = h->GetBinCenter(bm), s0 = 0.04 * keEl;
         TF1 g(TString::Format("gk%d_%d", i, b), "gaus", c0 - 2.5 * s0, c0 + 2.5 * s0);
         g.SetParameters(h->GetBinContent(bm), c0, s0);
         double mu = c0, sg = s0, emu = 0;
         if (h->Fit(&g, "QNR") == 0) {
            mu = g.GetParameter(1);
            sg = std::fabs(g.GetParameter(2));
            emu = g.GetParError(1);
         }
         double bias = mu - keEl;
         gKE[i].SetPoint(np[i], ctr, bias);
         gKE[i].SetPointError(np[i], 0, std::max(emu, sg / std::sqrt(N)));
         gRel[i].SetPoint(np[i], ctr, 100.0 * bias / keEl);
         gRel[i].SetPointError(np[i], 0, 100.0 * std::max(emu, sg / std::sqrt(N)) / keEl);
         gVsKE[i].SetPoint(np[i], keEl, 100.0 * bias / keEl);
         gVsKE[i].SetPointError(np[i], 0, 100.0 * std::max(emu, sg / std::sqrt(N)) / keEl);
         ++np[i];
         printf("  %8.3f %+7.3f %+6.2f%% %5.0f |", mu, bias, 100.0 * bias / keEl, N);
         delete h;
      }
      printf("\n");
   }

   TCanvas *c1 = new TCanvas("c1", "KE bias", 1500, 500);
   c1->Divide(3, 1);
   auto sty = [&](TGraphErrors &g, int i) {
      g.SetMarkerStyle(i == 0 ? 20 : 21);
      g.SetMarkerColor(col[i]);
      g.SetLineColor(col[i]);
      g.SetLineWidth(2);
      g.SetMarkerSize(1.2);
   };
   for (int i = 0; i < NS; ++i) {
      sty(gKE[i], i);
      sty(gRel[i], i);
      sty(gVsKE[i], i);
   }
   auto draw2 = [&](int pad, TGraphErrors *g, const char *title, double ylo, double yhi, double xlo, double xhi) {
      c1->cd(pad);
      g[0].SetTitle(title);
      g[0].GetYaxis()->SetRangeUser(ylo, yhi);
      g[0].GetXaxis()->SetLimits(xlo, xhi);
      g[0].Draw("ALP");
      g[1].Draw("LP same");
      auto *z = new TLine(xlo, 0, xhi, 0);
      z->SetLineStyle(2);
      z->SetLineColor(kGray + 2);
      z->Draw();
      auto *lg = new TLegend(0.62, 0.74, 0.89, 0.89);
      lg->AddEntry(&g[0], lbl[0], "lp");
      lg->AddEntry(&g[1], lbl[1], "lp");
      lg->Draw();
   };
   draw2(1, gKE, "KE bias vs #theta_{lab};#theta_{lab} [deg];KE_{fit} - KE_{elastic} [MeV]", -2.0, 1.0, thLo, thHi);
   draw2(2, gRel, "relative KE bias vs #theta_{lab};#theta_{lab} [deg];bias [%]", -12, 6, thLo, thHi);
   draw2(3, gVsKE, "relative KE bias vs elastic KE;KE_{elastic} [MeV];bias [%]", -12, 6, 0, 40);

   TString png = here + "/plots/ke_bias_C14.png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/ke_bias_C14.root", "RECREATE");
   for (int i = 0; i < NS; ++i) {
      gKE[i].Write(TString::Format("bias_%s", lbl[i]));
      gRel[i].Write(TString::Format("relbias_%s", lbl[i]));
      gVsKE[i].Write(TString::Format("relbias_vsKE_%s", lbl[i]));
   }
   fo.Close();
   printf("\nwrote %s\n\n", png.Data());
}
