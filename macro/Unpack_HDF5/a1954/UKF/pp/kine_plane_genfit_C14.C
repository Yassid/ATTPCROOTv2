/// @file kine_plane_genfit_C14.C
/// @brief Proton KE vs theta_lab, GENFIT only, with the two-body kinematic lines of 14C(p,p').
///
/// Single-panel version of kine_plane_C14.C for the write-up (Fig. 6.1): the adopted GENFIT+CATIMA
/// production, the adopted vertex slab, and the lines at the ANCHORED beam energy (159.75 MeV, from
/// the 6.094 MeV level -- not the 161.0 the older plot used). The theta is the measured one: the
/// empirical angle correction of the calibration chapter is NOT applied, so the plane shows what the
/// fitter delivers.
///
///   root -l -b -q 'kine_plane_genfit_C14.C()'

static double kg_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// Same invariant-mass expression as ex_C14.C (Eq. 6.3 of the write-up)
static double kg_ex(double m1, double m2, double Eb, double thl, double Ke)
{
   double m3 = m2, m4 = m1;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * (Eb + m1);
   double u = m2 * m2 + m3 * m3 - 2 * m2 * (Ke + m3);
   double a = (std::cos(thl) * kg_om2(s, m1 * m1, m2 * m2) * kg_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

/// All KE solutions giving this Ex at this lab angle (two branches near the kinematic ceiling)
static std::vector<double> kg_keFor(double m1, double m2, double Eb, double thl, double exWanted)
{
   std::vector<double> sol;
   auto g = [&](double k) { return kg_ex(m1, m2, Eb, thl, k) - exWanted; };
   const int NS = 4000;
   double k0 = 0.005, dk = (0.6 * Eb - k0) / NS;
   double f0 = g(k0);
   for (int i = 1; i <= NS; ++i) {
      double k1 = k0 + dk, f1 = g(k1);
      if (std::isfinite(f0) && std::isfinite(f1) && f0 * f1 <= 0) {
         double lo = k0, hi = k1, flo = f0;
         for (int j = 0; j < 60; ++j) {
            double mid = 0.5 * (lo + hi), fm = g(mid);
            if (flo * fm <= 0)
               hi = mid;
            else {
               lo = mid;
               flo = fm;
            }
         }
         sol.push_back(0.5 * (lo + hi));
      }
      k0 = k1;
      f0 = f1;
   }
   return sol;
}

void kine_plane_genfit_C14(TString file = "plots/proton_kin_cat5.root", Double_t Ebeam = 159.75,
                           TString cut = "vertexz>10&&vertexz<490", TString outPng = "plots/kine_plane_genfit_C14.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   TFile *f = TFile::Open(here + "/" + file);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mmissing %s\033[0m\n", file.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   auto *h = new TH2D("hk", "", 170, 5, 90, 200, 0, 40);
   Long64_t n = t->Draw("ke:theta>>hk", cut, "goff");
   h->SetDirectory(nullptr);
   printf("%s  cut \"%s\"  -> %lld of %lld tracks\n", file.Data(), cut.Data(), n, t->GetEntries());

   // inelastic loci have two KE solutions below their maximum lab angle: drawn as two graphs
   struct Level {
      double ex;
      const char *lbl;
      int col, ls;
   };
   std::vector<Level> levels = {{0.0, "g.s. 0^{+}", kRed + 1, 1},
                                {6.094, "6.094 1^{-}", kBlack, 1},
                                {8.317, "8.317 2^{+}", kMagenta + 2, 1}};

   TCanvas *c1 = new TCanvas("c1", "kinematics plane GENFIT", 1600, 1240);
   c1->SetLogz();
   c1->SetRightMargin(0.13);
   c1->SetLeftMargin(0.11);
   h->SetTitle(";#theta_{lab} [deg];K_{p} [MeV]");
   h->GetYaxis()->SetTitleOffset(1.1);
   h->Draw("colz");

   auto *lg = new TLegend(0.52, 0.66, 0.85, 0.88);
   lg->SetHeader(TString::Format("^{14}C(p,p'), E_{beam} = %.2f MeV", Ebeam));
   lg->SetTextSize(0.034);
   for (auto &L : levels) {
      auto *lo = new TGraph(), *hi = new TGraph();
      for (double th = 5; th <= 89.9; th += 0.1) {
         auto s = kg_keFor(m1, m2, Ebeam, th * TMath::DegToRad(), L.ex);
         if (s.empty())
            continue;
         std::sort(s.begin(), s.end());
         hi->SetPoint(hi->GetN(), th, s.back());
         if (s.size() > 1)
            lo->SetPoint(lo->GetN(), th, s.front());
      }
      for (auto *g : {hi, lo}) {
         g->SetLineColor(L.col);
         g->SetLineWidth(5);
         g->SetLineStyle(L.ls);
         if (g->GetN() > 1)
            g->Draw("L same");
      }
      lg->AddEntry(hi, TString::Format("E_{x} = %s", L.lbl), "l");
      printf("  Ex %6.3f : %d pts upper branch, %d lower; theta_max %.1f deg\n", L.ex, hi->GetN(), lo->GetN(),
             hi->GetN() ? hi->GetPointX(hi->GetN() - 1) : NAN);
   }
   lg->Draw();

   TString png = here + "/" + outPng;
   c1->SaveAs(png);
   c1->SaveAs(TString(png).ReplaceAll(".png", ".pdf"));
   printf("wrote %s\n", png.Data());
}
