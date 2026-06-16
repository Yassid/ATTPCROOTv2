/// @file ex_kin_brho_dp.C
/// @brief Proton KE-from-Brho vs theta_lab kinematics for 16C(d,p)17C. Same as
///        ex_kin_dp.C but the energy axis is derived from the Spyral/PRA magnetic
///        rigidity Brho (an INDEPENDENT momentum estimate from the PID plane, not the
///        genfit fitted KE) -> a cross-check that the two momentum measurements agree.
///        For a proton (Z=1): p[MeV/c] = 299.792458 * Brho[T m];
///        KE = sqrt(p^2 + m_p^2) - m_p.  Reads the cached proton_kin_dp.root ntuple.
///
///   root -l -b -q 'ex_kin_brho_dp.C'            // all angles
///   root -l -b -q 'ex_kin_brho_dp.C(192.0,50.0,"proton_kin_dp.root",90,180)'  // backward only

#include <vector>

static double omega2_b(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double kine_ex_b(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * omega2_b(s, m1 * m1, m2 * m2) * omega2_b(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}
static TGraph *kin_curve_b(double m1, double m2, double m3, double m4, double Ebeam, double exTarget, double thLoDeg,
                           double thHiDeg, double keMax)
{
   TGraph *g = new TGraph();
   int np = 0;
   for (double thd = thLoDeg; thd <= thHiDeg; thd += 0.5) {
      double th = thd * TMath::DegToRad();
      double prevE = kine_ex_b(m1, m2, m3, m4, Ebeam, th, 0.02), prevK = 0.02;
      for (double ke = 0.04; ke <= keMax; ke += 0.02) {
         double e = kine_ex_b(m1, m2, m3, m4, Ebeam, th, ke);
         if (prevE > -1e8 && e > -1e8 && (prevE - exTarget) * (e - exTarget) <= 0)
            g->SetPoint(np++, thd, prevK + (exTarget - prevE) / (e - prevE) * (ke - prevK));
         prevE = e;
         prevK = ke;
      }
   }
   return g;
}

void ex_kin_brho_dp(double Ebeam = 192.0, double keMax = 50.0, const char *cache = "proton_kin_dp.root",
                    double thLoCut = 0.0, double thHiCut = 180.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.01410178 * u, m_p = 1.00782503 * u, m_C17 = 17.0225864 * u;

   TFile *f = TFile::Open(cache);
   TNtuple *t = (TNtuple *)f->Get("pk");
   if (!t) {
      printf("no ntuple 'pk' in %s\n", cache);
      return;
   }

   TH2F *hk = new TH2F("hkb", "^{16}C(d,p)^{17}C kinematics from B#rho;#theta_{lab} [deg];proton KE(B#rho) [MeV]", 180,
                       0, 180, 150, 0, 25);
   // KE from Brho on the fly: p = 299.792458*brho (MeV/c, Z=1); KE = sqrt(p^2+mp^2)-mp
   TString ke_brho = "(sqrt(pow(299.792458*brho,2)+938.272*938.272)-938.272)";
   t->Draw(ke_brho + ":theta>>hkb", Form("theta>=%g&&theta<=%g", thLoCut, thHiCut), "goff");

   double exVals[] = {0.0, 5.0, 10.0};
   Color_t cols[] = {kRed, kAzure + 1, kGreen + 2};
   TGraph *gc[3];
   for (int i = 0; i < 3; ++i) {
      gc[i] = kin_curve_b(m_C16, m_d, m_p, m_C17, Ebeam, exVals[i], 1, 170, keMax);
      gc[i]->SetLineColor(cols[i]);
      gc[i]->SetMarkerColor(cols[i]);
      gc[i]->SetMarkerStyle(1);
   }

   TCanvas *c = new TCanvas("cb", "kin_brho", 1100, 800);
   c->SetRightMargin(0.13);
   hk->Draw("colz");
   for (int i = 0; i < 3; ++i)
      gc[i]->Draw("P same");
   TLegend *leg = new TLegend(0.6, 0.70, 0.86, 0.88);
   leg->SetHeader(Form("two-body, E_{beam}=%.0f MeV", Ebeam));
   for (int i = 0; i < 3; ++i)
      leg->AddEntry(gc[i], Form("E_{x} = %.0f MeV", exVals[i]), "p");
   leg->Draw();
   c->SaveAs("plots/kin_brho_dp.png");
   printf("saved plots/kin_brho_dp.png  (energy axis = KE derived from Brho, not genfit KE)\n");
}
