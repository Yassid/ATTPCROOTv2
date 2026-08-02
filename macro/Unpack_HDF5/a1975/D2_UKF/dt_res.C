/// @file dt_res.C
/// @brief Where does the 16C(d,t)15C excitation-energy resolution actually come from?
///
/// The integrated (d,t) spectrum is far broader than the 16C(p,d)15C one (sigma 0.25 MeV,
/// FWHM 0.59, from a1975_panels/pd_tight_panels.root) even though both channels reconstruct
/// the SAME residual nucleus from the same beam. This macro takes the ex_dt_a1975.C cache
/// apart along the three axes that can carry the difference:
///
///   1. theta_lab. The (d,t) g.s. locus turns over at theta_max ~ 56 deg, and the old
///      2 < KE_t < 20 MeV box sits right on that turnover. There dKE/dtheta -> infinity,
///      so a fixed angular error maps onto an ever-larger Ex error: the resolution is
///      expected to DEGRADE towards the turnover and be best on the forward, high-KE part
///      of the locus that the KE box throws away. sigma(theta) tests exactly that.
///   2. vertex z. An active target reconstructs the missing mass with ONE beam energy for
///      every event, but the 16C loses energy continuously in the gas, so Ex acquires a
///      slope in z. The (p,d) control measures +0.022 MeV/100 mm (0.18 MeV over the whole
///      drift), small next to its own sigma -- if (d,t) is much steeper, that is the cause.
///   3. Ebeam. A wrong beam energy is a SCALE error, so it tilts Ex against theta rather
///      than shifting it. The tilt is what broadens the theta-integrated peak. The (p,d)
///      control sits at -0.056 MeV/10 deg with Ebeam = 195.5; the (d,t) cache is built at
///      Ebeam = 180 (inherited from the old C16_dt_anaFit.C, where 192 was commented out).
///
/// Ex is RECOMPUTED here from (ke, theta) so Ebeam can be changed without rebuilding the
/// cache. Pass ebeamScan>0 to sweep it and watch the theta tilt cross zero.
///
///   root -b -q 'dt_res.C("/path/dt_kin_full.root")'
///   root -b -q 'dt_res.C("/path/dt_kin_full.root", 195.5)'

static double om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// missing-mass Ex of the residual, identical form to ex_dt_a1975.C's kine_2b
static double exOf(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, u = m2 * m2 + m3 * m3 - 2 * m2 * E3;
   double a = (std::cos(th) * om2(s, m1 * m1, m2 * m2) * om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a < 0 ? std::nan("") : std::sqrt(a) - m4;
}

/// gaussian + linear background around a seeded centroid; returns false if the fit is not usable
static bool fitPk(TH1F *h, double c, double win, double &mu, double &er, double &sg)
{
   if (h->GetEntries() < 60)
      return false;
   TF1 g("g", "gaus(0)+pol1(3)", c - win, c + win);
   int bm = h->GetMaximumBin();
   g.SetParameters(h->GetBinContent(bm), c, 0.35, 0, 0);
   g.SetParLimits(1, c - win, c + win);
   g.SetParLimits(2, 0.05, 2.0);
   if (h->Fit(&g, "QRN") != 0)
      return false;
   mu = g.GetParameter(1);
   er = g.GetParError(1);
   sg = std::fabs(g.GetParameter(2));
   return er > 0 && er < 1.0;
}

void dt_res(TString cache = "/tmp/dt_kin_full.root", double Ebeam = 180.0, double icLo = 900, double icHi = 1300,
            double chi2max = 5.0, double keLo = 2.0, double keHi = 90.0, double ebeamScan = 0,
            TString plotOut = "plots/dt_res.png", double thLo = 0.0, double thHi = 180.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double u = 931.49401;
   const double m1 = 16.0147013 * u, m2 = 2.0135532 * u, m3 = 3.01550072 * u, m4 = 15.0105993 * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) {
      printf("no tree `pk` in %s\n", cache.Data());
      return;
   }
   float ke, theta, vertexz, chi2ndf, ic = -1, thcm = 0, run = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &theta);
   t->SetBranchAddress("vertexz", &vertexz);
   t->SetBranchAddress("chi2ndf", &chi2ndf);
   if (t->GetBranch("thcm"))
      t->SetBranchAddress("thcm", &thcm);
   if (t->GetBranch("run"))
      t->SetBranchAddress("run", &run);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC)
      t->SetBranchAddress("ic", &ic);

   struct Ev {
      float ke, th, vz, ic, run;
   };
   std::vector<Ev> ev;
   long nAll = 0, nChi = 0, nIC = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      ++nAll;
      if (chi2ndf > chi2max || ke <= keLo || ke > keHi)
         continue;
      if (theta < thLo || theta >= thHi)
         continue;
      ++nChi;
      if (hasIC && icLo > 0 && (ic < icLo || ic > icHi))
         continue;
      ++nIC;
      ev.push_back({ke, theta, vertexz, ic, run});
   }
   printf("\n=== dt_res: %s ===\n", gSystem->BaseName(cache));
   printf("cache %ld -> chi2/KE %ld -> IC[%.0f,%.0f] %ld candidates   (Ebeam %.1f)\n", nAll, nChi, icLo, icHi, nIC,
          Ebeam);
   if (ev.empty()) {
      printf("nothing survives the cuts\n");
      return;
   }

   auto exAt = [&](const Ev &e, double E) { return exOf(m1, m2, m3, m4, E, e.th * TMath::DegToRad(), e.ke); };

   // ---------------------------------------------------------------- 1. sigma vs theta
   printf("\n-- resolution vs theta_lab (g.s. region, Ebeam %.1f) --\n", Ebeam);
   printf("%-14s %8s %9s %9s %9s\n", "theta window", "N", "mu", "sigma", "FWHM");
   const double thEdge[] = {20, 25, 30, 35, 40, 45, 50, 56, 65, 90};
   const int nTh = sizeof(thEdge) / sizeof(double) - 1;
   TGraph *gSig = new TGraph();
   for (int i = 0; i < nTh; ++i) {
      TH1F h(Form("hth%d", i), "", 240, -8, 16);
      h.SetDirectory(nullptr);
      for (auto &e : ev)
         if (e.th >= thEdge[i] && e.th < thEdge[i + 1]) {
            double x = exAt(e, Ebeam);
            if (!std::isnan(x))
               h.Fill(x);
         }
      // seed at the bin's OWN maximum: the theta tilt moves the peak far from Ex=0, and a
      // fit window pinned at zero then latches onto whatever shoulder is nearest instead.
      double seed = h.GetEntries() > 0 ? h.GetBinCenter(h.GetMaximumBin()) : 0.0;
      double mu, er, sg;
      bool ok = fitPk(&h, seed, 1.2, mu, er, sg);
      printf("%5.0f - %-6.0f %8.0f %9s %9s %9s\n", thEdge[i], thEdge[i + 1], h.GetEntries(),
             ok ? Form("%+.3f", mu) : "-", ok ? Form("%.3f", sg) : "-", ok ? Form("%.3f", 2.355 * sg) : "-");
      if (ok)
         gSig->SetPoint(gSig->GetN(), 0.5 * (thEdge[i] + thEdge[i + 1]), sg);
   }

   // ---------------------------------------------------------------- 2. Ex vs vertex z
   auto *pvz = new TProfile("pvz", "E_{x} vs vertex z;vertex z [mm];#LTE_{x}#GT [MeV]", 20, 0, 1000);
   auto *pth = new TProfile("pth", "E_{x} vs #theta_{lab};#theta_{lab} [deg];#LTE_{x}#GT [MeV]", 22, 15, 70);
   for (auto &e : ev) {
      double x = exAt(e, Ebeam);
      if (std::isnan(x) || std::fabs(x) > 3)
         continue; // g.s. region only
      pvz->Fill(e.vz, x);
      pth->Fill(e.th, x);
   }
   TF1 lz("lz", "pol1", 100, 900);
   pvz->Fit(&lz, "QRN");
   TF1 lt("lt", "pol1", 25, 60);
   pth->Fit(&lt, "QRN");
   printf("\n-- drifts in the g.s. region (|Ex|<3) --\n");
   printf("Ex vs vertex z : %+.4f MeV / 100 mm  (%.2f MeV across 800 mm)\n", 100 * lz.GetParameter(1),
          800 * lz.GetParameter(1));
   printf("Ex vs theta    : %+.4f MeV / 10 deg   (%.2f MeV across 25-60 deg)\n", 10 * lt.GetParameter(1),
          35 * lt.GetParameter(1));
   printf("   (p,d) control: +0.0220 MeV/100 mm, -0.0559 MeV/10 deg, sigma 0.252)\n");

   // ---------------------------------------------------------------- 3. optional Ebeam sweep
   if (ebeamScan > 0) {
      printf("\n-- Ebeam sweep: theta tilt and g.s. width --\n");
      printf("%-9s %11s %9s %9s %8s\n", "Ebeam", "tilt/10deg", "gs_mu", "gs_sig", "N");
      for (double E = Ebeam - ebeamScan; E <= Ebeam + ebeamScan + 1e-9; E += ebeamScan / 5) {
         TProfile p("p", "", 22, 15, 70);
         p.SetDirectory(nullptr);
         TH1F h("h", "", 240, -8, 16);
         h.SetDirectory(nullptr);
         for (auto &e : ev) {
            double x = exAt(e, E);
            if (std::isnan(x))
               continue;
            h.Fill(x);
            if (std::fabs(x) < 3)
               p.Fill(e.th, x);
         }
         TF1 l("l", "pol1", 25, 60);
         p.Fit(&l, "QRN");
         double mu, er, sg;
         bool ok = fitPk(&h, 0.0, 2.0, mu, er, sg);
         printf("%-9.1f %11.4f %9s %9s %8.0f\n", E, 10 * l.GetParameter(1), ok ? Form("%+.3f", mu) : "-",
                ok ? Form("%.3f", sg) : "-", h.GetEntries());
      }
   }

   // ---------------------------------------------------------------- plots
   auto *hExTh = new TH2F("hExTh", "E_{x} vs #theta_{lab};#theta_{lab} [deg];E_{x}(^{15}C) [MeV]", 70, 15, 85, 160,
                          -8, 16);
   auto *hKeTh = new TH2F("hKeTh", "KE_{t} vs #theta_{lab} (gated);#theta_{lab} [deg];KE_{t} [MeV]", 80, 10, 90, 180,
                          0, 90);
   for (auto &e : ev) {
      double x = exAt(e, Ebeam);
      if (!std::isnan(x))
         hExTh->Fill(e.th, x);
      hKeTh->Fill(e.th, e.ke);
   }
   auto *cv = new TCanvas("cvres", "dt_res", 1700, 1000);
   cv->Divide(3, 2);
   cv->cd(1);
   gPad->SetLogz();
   hKeTh->Draw("colz");
   cv->cd(2);
   gPad->SetLogz();
   hExTh->Draw("colz");
   for (double lv : {0.0, 0.740, 3.103}) {
      auto *l = new TLine(15, lv, 85, lv);
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
   }
   cv->cd(3);
   gSig->SetMarkerStyle(20);
   gSig->SetTitle("g.s. #sigma vs #theta_{lab};#theta_{lab} [deg];#sigma(E_{x}) [MeV]");
   gSig->Draw("APL");
   {
      auto *r = new TLine(20, 0.252, 90, 0.252);
      r->SetLineColor(kBlue);
      r->SetLineStyle(2);
      r->Draw();
   }
   cv->cd(4);
   pvz->Draw();
   cv->cd(5);
   pth->Draw();
   cv->cd(6);
   auto *hAll = new TH1F("hAll", "E_{x}, all gated;E_{x}(^{15}C) [MeV];tritons", 160, -8, 16);
   for (auto &e : ev) {
      double x = exAt(e, Ebeam);
      if (!std::isnan(x))
         hAll->Fill(x);
   }
   hAll->Draw();
   for (double lv : {0.0, 0.740, 3.103, 4.220}) {
      auto *l = new TLine(lv, 0, lv, hAll->GetMaximum());
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
   }
   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("\nsaved %s\n", plotOut.Data());
}
