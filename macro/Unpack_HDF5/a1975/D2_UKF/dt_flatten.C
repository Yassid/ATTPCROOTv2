/// @file dt_flatten.C
/// @brief The 16C(d,t)15C peak is sharp in every theta slice and smeared only after the
///        slices are added up. Find what makes the slices disagree.
///
/// dt_res.C measured, on the slow branch (45 < theta_lab < 56, the region the old
/// C16_dt_anaFit.C used): sigma = 0.21 MeV inside 45-50 deg, but a centroid that moves from
/// +3.37 to +0.52 MeV between 45-50 and 50-56. Per slice the resolution already matches the
/// 16C(p,d)15C control (sigma 0.252); it is the ~1.7 MeV/10 deg drift ACROSS slices that
/// turns the summed spectrum into a blob. So the question is not "why is genfit worse here"
/// but "what scale error tilts Ex against theta".
///
/// Two candidates produce exactly this tilt, and they are separable because they tilt it
/// differently:
///   Ebeam  -- the missing mass is reconstructed with one beam energy for every event; too
///             low a value makes Ex fall with theta (the (d,t) cache is built at 180 MeV,
///             inherited from C16_dt_anaFit.C where 192 sits commented out, while the
///             (p,p') 2+ ruler put the H2-run beam at 195.5).
///   keScale -- a multiplicative error in the fitted triton energy, e.g. from the missing
///             energy-loss correction (the D2 genfit production runs matEffects = kFALSE
///             and passes no dE/dx table, so nothing puts back what the triton lost in the
///             gas between the vertex and the pad plane).
///
/// The metric is deliberately fit-free, because a background-dominated slice defeats any
/// seeded gaussian: build the Ex histogram, take the tallest bin after light smoothing, and
/// measure its FWHM by walking down both sides. A tilt that is removed shows up as the
/// summed peak getting taller and narrower at the same time.
///
///   root -b -q 'dt_flatten.C("/path/dt_kin_full.root")'

static double om2f(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double exOfF(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, u = m2 * m2 + m3 * m3 - 2 * m2 * E3;
   double a = (std::cos(th) * om2f(s, m1 * m1, m2 * m2) * om2f(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a < 0 ? std::nan("") : std::sqrt(a) - m4;
}

/// tallest-structure width, no seed and no fit: returns FWHM and fills peak position/height
static double fwhmOf(TH1F *h, double &pos, double &hgt)
{
   TH1F s(*h);
   s.SetDirectory(nullptr);
   s.Smooth(2);
   int bm = s.GetMaximumBin();
   double half = 0.5 * s.GetBinContent(bm);
   pos = s.GetBinCenter(bm);
   hgt = s.GetBinContent(bm);
   int lo = bm, hi = bm;
   while (lo > 1 && s.GetBinContent(lo) > half)
      --lo;
   while (hi < s.GetNbinsX() && s.GetBinContent(hi) > half)
      ++hi;
   return s.GetBinCenter(hi) - s.GetBinCenter(lo);
}

void dt_flatten(TString cache = "/tmp/dt_kin_full.root", double icLo = 900, double icHi = 1300, double chi2max = 5.0,
                double thLo = 45, double thHi = 56, double keLo = 2, double keHi = 20,
                double eLo = 160, double eHi = 215, double dE = 5, double ksLo = 0.90, double ksHi = 1.30,
                double dks = 0.05, TString plotOut = "plots/dt_flatten.png")
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = 16.0147013 * u, m2 = 2.0135532 * u, m3 = 3.01550072 * u, m4 = 15.0105993 * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   float ke, theta, vertexz, chi2ndf, ic = -1;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &theta);
   t->SetBranchAddress("vertexz", &vertexz);
   t->SetBranchAddress("chi2ndf", &chi2ndf);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC)
      t->SetBranchAddress("ic", &ic);

   struct Ev {
      float ke, th, vz;
   };
   std::vector<Ev> ev;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (chi2ndf > chi2max || ke <= keLo || ke > keHi)
         continue;
      if (theta < thLo || theta >= thHi)
         continue;
      if (hasIC && icLo > 0 && (ic < icLo || ic > icHi))
         continue;
      ev.push_back({ke, theta, vertexz});
   }
   printf("\n=== dt_flatten: %zu events, theta[%.0f,%.0f] KE[%.0f,%.0f] ===\n", ev.size(), thLo, thHi, keLo, keHi);
   if (ev.size() < 100) {
      printf("too few events\n");
      return;
   }

   auto build = [&](double E, double ks) {
      TH1F *h = new TH1F("hb", "", 200, -6, 14);
      h->SetDirectory(nullptr);
      for (auto &e : ev) {
         double x = exOfF(m1, m2, m3, m4, E, e.th * TMath::DegToRad(), ks * e.ke);
         if (!std::isnan(x))
            h->Fill(x);
      }
      return h;
   };

   // ---- 1D: Ebeam alone --------------------------------------------------------------
   printf("\n-- Ebeam alone (keScale = 1) --\n");
   printf("%-8s %9s %9s %9s\n", "Ebeam", "peak", "FWHM", "height");
   double bestE = eLo, bestW = 1e9;
   TGraph *gE = new TGraph();
   for (double E = eLo; E <= eHi + 1e-9; E += dE) {
      TH1F *h = build(E, 1.0);
      double p, ht, w = fwhmOf(h, p, ht);
      printf("%-8.1f %9.3f %9.3f %9.0f\n", E, p, w, ht);
      gE->SetPoint(gE->GetN(), E, w);
      if (w < bestW) {
         bestW = w;
         bestE = E;
      }
      delete h;
   }
   printf("  narrowest at Ebeam = %.1f  (FWHM %.3f)\n", bestE, bestW);

   // ---- 2D: Ebeam x keScale ----------------------------------------------------------
   printf("\n-- Ebeam x keScale grid, FWHM of the tallest structure --\n");
   printf("%-8s", "Ebeam\\ks");
   for (double k = ksLo; k <= ksHi + 1e-9; k += dks)
      printf(" %6.2f", k);
   printf("\n");
   double bE = 0, bK = 0, bW = 1e9;
   auto *h2 = new TH2F("h2", "FWHM of the summed peak;E_{beam} [MeV];KE scale", int((eHi - eLo) / dE) + 1,
                       eLo - dE / 2, eHi + dE / 2, int((ksHi - ksLo) / dks) + 1, ksLo - dks / 2, ksHi + dks / 2);
   for (double E = eLo; E <= eHi + 1e-9; E += dE) {
      printf("%-8.1f", E);
      for (double k = ksLo; k <= ksHi + 1e-9; k += dks) {
         TH1F *h = build(E, k);
         double p, ht, w = fwhmOf(h, p, ht);
         printf(" %6.2f", w);
         h2->Fill(E, k, w);
         if (w < bW) {
            bW = w;
            bE = E;
            bK = k;
         }
         delete h;
      }
      printf("\n");
   }
   printf("  best: Ebeam %.1f  keScale %.2f  FWHM %.3f\n", bE, bK, bW);

   // ---- per-theta centroids, before and after ----------------------------------------
   auto muVsTheta = [&](double E, double ks, const char *tag) {
      printf("\n-- centroid vs theta at Ebeam %.1f, keScale %.2f (%s) --\n", E, ks, tag);
      double lo = 1e9, hi = -1e9;
      for (double a = thLo; a < thHi - 1e-9; a += 2.0) {
         TH1F h("hs", "", 200, -6, 14);
         h.SetDirectory(nullptr);
         for (auto &e : ev)
            if (e.th >= a && e.th < a + 2.0) {
               double x = exOfF(m1, m2, m3, m4, E, e.th * TMath::DegToRad(), ks * e.ke);
               if (!std::isnan(x))
                  h.Fill(x);
            }
         if (h.GetEntries() < 40) {
            printf("  %4.0f-%-4.0f  %5.0f      (too few)\n", a, a + 2, h.GetEntries());
            continue;
         }
         double p, ht, w = fwhmOf(&h, p, ht);
         printf("  %4.0f-%-4.0f  %5.0f   peak %+6.3f   FWHM %5.3f\n", a, a + 2, h.GetEntries(), p, w);
         lo = std::min(lo, p);
         hi = std::max(hi, p);
      }
      if (hi > lo)
         printf("  spread of the slice centroids: %.3f MeV\n", hi - lo);
   };
   muVsTheta(180.0, 1.0, "as-produced");
   muVsTheta(bE, bK, "best grid point");

   auto *cv = new TCanvas("cvf", "dt_flatten", 1500, 500);
   cv->Divide(3, 1);
   cv->cd(1);
   gE->SetMarkerStyle(20);
   gE->SetTitle("summed-peak FWHM vs E_{beam};E_{beam} [MeV];FWHM [MeV]");
   gE->Draw("APL");
   cv->cd(2);
   h2->Draw("colz");
   cv->cd(3);
   TH1F *hb = build(bE, bK);
   hb->SetTitle(Form("E_{x} at E_{beam}=%.1f, keScale=%.2f;E_{x}(^{15}C) [MeV];tritons", bE, bK));
   hb->Draw();
   TH1F *h0 = build(180.0, 1.0);
   h0->SetLineColor(kGray + 2);
   h0->Draw("same");
   for (double lv : {0.0, 0.740, 3.103, 4.220}) {
      auto *l = new TLine(lv, 0, lv, hb->GetMaximum());
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
   }
   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("\nsaved %s\n", plotOut.Data());
}
