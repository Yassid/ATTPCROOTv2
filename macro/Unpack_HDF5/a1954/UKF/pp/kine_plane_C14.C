/// @file kine_plane_C14.C
/// @brief Proton KE vs theta_lab with the elastic and inelastic kinematic lines drawn on top.
///
/// This is the plot of Fig. 5 (upper) of Ayyad et al., EPJ A 59:294 (2023), for the same data.
/// Everything else in this study is a projection of it: the E_x locus drift, the KE bias, the
/// missing FRESCO secondary maximum. Seeing where the measured ridge leaves the kinematic line,
/// and in which direction, says what kind of error it is:
///   * ridge parallel to the line but offset -> beam energy
///   * ridge crossing the line with a tilt   -> angle or momentum scale
///   * ridge bending away only at high KE    -> energy loss / curvature at large radius
///
/// Lines are computed by bisecting the same two-body expression the analysis uses, at the given
/// beam energy, for E_x = 0 and for the 6.094 MeV state.
///
///   root -b -q 'kine_plane_C14.C()'

#include <tuple>

static double kp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double kp_ex(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * kp_om2(s, m1 * m1, m2 * m2) * kp_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

/// KE giving this Ex at this lab angle
static double kp_keFor(double m1, double m2, double Eb, double thl, double exWanted)
{
   double lo = 0.02, hi = 0.999 * Eb;
   auto g = [&](double k) { return kp_ex(m1, m2, m2, m1, Eb, thl, k) - exWanted; };
   double flo = g(lo), fhi = g(hi);
   if (!std::isfinite(flo) || !std::isfinite(fhi) || flo * fhi > 0)
      return NAN;
   for (int i = 0; i < 120; ++i) {
      double mid = 0.5 * (lo + hi), fm = g(mid);
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

void kine_plane_C14(Double_t Eb1 = 161.0, Double_t Eb2 = 168.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   const int NS = 2;
   const char *file[NS] = {"plots/proton_kin_300_ukf_nc.root", "plots/proton_kin_300gfx_nc.root"};
   const char *lbl[NS] = {"UKF", "GENFIT"};

   // kinematic lines
   auto line = [&](double Eb, double ex, int col, int ls) {
      auto *g = new TGraph();
      int n = 0;
      for (double th = 5; th <= 89.5; th += 0.5) {
         double ke = kp_keFor(m1, m2, Eb, th * TMath::DegToRad(), ex);
         if (std::isfinite(ke) && ke > 0)
            g->SetPoint(n++, th, ke);
      }
      g->SetLineColor(col);
      g->SetLineWidth(3);
      g->SetLineStyle(ls);
      return g;
   };
   TGraph *el1 = line(Eb1, 0.0, kRed + 1, 1);
   TGraph *el2 = line(Eb2, 0.0, kMagenta + 2, 2);
   TGraph *in1 = line(Eb1, 6.094, kBlack, 1);
   TGraph *in2 = line(Eb1, 7.012, kGray + 2, 7);

   TCanvas *c1 = new TCanvas("c1", "kinematics plane", 1500, 1000);
   c1->Divide(2, 2);

   for (int i = 0; i < NS; ++i) {
      TFile *f = TFile::Open(here + "/" + file[i]);
      if (!f || f->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", file[i]);
         return;
      }
      TTree *t = (TTree *)f->Get("pk");
      auto *h = new TH2D(TString::Format("hk%d", i), "", 170, 5, 90, 200, 0, 40);
      t->Draw(TString::Format("ke:theta>>hk%d", i), "", "goff");
      h->SetDirectory(nullptr);

      c1->cd(i + 1);
      gPad->SetLogz();
      h->SetTitle(TString::Format("%s: measured kinematics;#theta_{lab} [deg];KE [MeV]", lbl[i]));
      h->Draw("colz");
      el1->Draw("L same");
      el2->Draw("L same");
      in1->Draw("L same");
      in2->Draw("L same");
      if (i == 0) {
         auto *lg = new TLegend(0.50, 0.62, 0.89, 0.88);
         lg->AddEntry(el1, TString::Format("elastic, E_{beam}=%.0f", Eb1), "l");
         lg->AddEntry(el2, TString::Format("elastic, E_{beam}=%.0f", Eb2), "l");
         lg->AddEntry(in1, "E_{x}=6.094, E_{beam}=161", "l");
         lg->AddEntry(in2, "E_{x}=7.012, E_{beam}=161", "l");
         lg->SetTextSize(0.030);
         lg->Draw();
      }

      // measured ridge: mode of KE in each theta slice, restricted to the elastic band
      auto *gr = new TGraph();
      int n = 0;
      printf("\n===== %s : elastic ridge vs the kinematic line =====\n", lbl[i]);
      printf("  theta_lab |  KE(Eb=%.0f)  KE(Eb=%.0f) |  measured ridge   ridge-line(%.0f)   rel\n", Eb1, Eb2, Eb1);
      for (double th = 22; th < 84; th += 3) {
         double keL1 = kp_keFor(m1, m2, Eb1, (th + 1.5) * TMath::DegToRad(), 0.0);
         double keL2 = kp_keFor(m1, m2, Eb2, (th + 1.5) * TMath::DegToRad(), 0.0);
         if (!std::isfinite(keL1))
            continue;
         auto *hs = new TH1D("hs", "", 160, std::max(0.0, 0.55 * keL1), 1.55 * keL1);
         t->Draw("ke>>hs", TString::Format("theta>=%g&&theta<%g", th, th + 3), "goff");
         hs->SetDirectory(nullptr);
         if (hs->Integral() > 100) {
            hs->Smooth(1);
            double pk = hs->GetBinCenter(hs->GetMaximumBin());
            gr->SetPoint(n++, th + 1.5, pk);
            printf("  %4.0f-%-4.0f | %9.2f %11.2f |  %12.2f %14.2f %8.1f%%\n", th, th + 3, keL1, keL2, pk, pk - keL1,
                   100 * (pk - keL1) / keL1);
         }
         delete hs;
      }
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kWhite);
      gr->SetMarkerSize(0.9);
      gr->SetLineColor(kWhite);
      gr->SetLineWidth(2);
      gr->Draw("LP same");

      c1->cd(i + 3);
      auto *hd = new TH2D(TString::Format("hd%d", i), "", 170, 5, 90, 160, -8, 8);
      // KE minus the elastic line, per track
      float ke, th;
      t->SetBranchAddress("ke", &ke);
      t->SetBranchAddress("theta", &th);
      // precompute the line on a fine grid for speed
      const int NG = 850;
      std::vector<double> grid(NG, NAN);
      for (int g = 0; g < NG; ++g)
         grid[g] = kp_keFor(m1, m2, Eb1, (5.0 + 0.1 * g) * TMath::DegToRad(), 0.0);
      Long64_t N = t->GetEntries();
      for (Long64_t e = 0; e < N; ++e) {
         t->GetEntry(e);
         int g = (int)((th - 5.0) / 0.1);
         if (g < 0 || g >= NG || !std::isfinite(grid[g]))
            continue;
         hd->Fill(th, ke - grid[g]);
      }
      hd->SetDirectory(nullptr);
      gPad->SetLogz();
      hd->SetTitle(TString::Format("%s: KE - KE_{elastic}(E_{beam}=%.0f);#theta_{lab} [deg];#DeltaKE [MeV]", lbl[i],
                                   Eb1));
      hd->Draw("colz");
      auto *z = new TLine(5, 0, 90, 0);
      z->SetLineStyle(2);
      z->SetLineColor(kWhite);
      z->SetLineWidth(2);
      z->Draw();
      f->Close();
   }

   TString png = here + "/plots/kine_plane_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
