/// @file dt_dvscan.C
/// @brief Find the a1975 D2 drift velocity from the data, on the one observable a two-body
///        reaction fixes exactly: E_x cannot depend on theta_cm.
///
/// The D2 chain runs at DriftVelocity = 1.15 cm/us (parameters/ATTPC.a1975_deuterium.par); the
/// H2 (p,d) benchmark runs at 1.30 (ATTPC.a1954.par). A wrong dv is a pure scale error on the
/// drift coordinate, so it can be swept on an EXISTING cache without re-running reco:
///
///   the pad plane fixes the transverse coordinates, so the fitted helix radius -- and with it
///   p_T = 0.3*B*R -- is untouched by dv. Only the pitch moves. With k = dv_true/dv_used,
///        p_z -> k*p_z ,  p_T invariant
///        theta' = atan2(p_T, k*p_z) ,  KE' = sqrt(p_T^2 + k^2 p_z^2 + m^2) - m
///   and the vertex z scales by k too. atan2 keeps backward tracks in their own hemisphere.
///
/// USE THE ELASTIC CHANNEL. 16C(d,d)16C has E_x = 0 identically, at every angle, with no level
/// scheme to argue about and no weak peak to fit -- exactly what a scale calibration wants. The
/// (d,t) spectrum is far too thin for this: fitting its 3.103 level while dv moves the whole
/// scale underneath makes the fit change which structure it has hold of, and the answer wanders
/// between parameter limits. So the default cache here is deuteron_kin_dd.root.
///
/// Everything is fit-free on purpose -- tallest bin of a lightly smoothed histogram, and FWHM by
/// walking down its sides. A seeded gaussian is what broke the first version of this scan.
///
/// Two conditions, one parameter, so the problem is over-determined and their agreement is a
/// real test rather than a fit:
///   (1) the elastic peak must be FLAT across theta_cm
///   (2) the elastic peak must sit AT ZERO
///
///   root -b -q 'dt_dvscan.C()'                                  // (d,d) elastic, the calibration
///   root -b -q 'dt_dvscan.C("/mnt/f/a1975/dt_kin_full.root",1.15,0.85,1.35,0.05,180,3.01550072,15.0105993,3.103)'
///
static double om2d(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static std::pair<double, double> kine2bd(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3, E4 = E1 + m2 - E3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, u = m2 * m2 + m3 * m3 - 2 * m2 * E3;
   double a = (std::cos(th) * om2d(s, m1 * m1, m2 * m2) * om2d(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   if (a < 0)
      return {std::nan(""), std::nan("")};
   double m4x = std::sqrt(a), Ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * E4;
   double arg = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
                 (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
                (om2d(s, m1 * m1, m2 * m2) * om2d(s, m3 * m3, m4x * m4x));
   if (arg < -1 || arg > 1)
      return {Ex, std::nan("")};
   return {Ex, (TMath::Pi() - std::acos(arg)) * TMath::RadToDeg()};
}

/// tallest bin of a lightly smoothed copy, and its FWHM by walking down both sides. No seed,
/// no parameter limits, nothing that can latch onto the wrong structure as the scale moves.
static double peakOf(TH1F *h, double &fwhm, double lo = -1e9, double hi = 1e9)
{
   TH1F s(*h);
   s.SetDirectory(nullptr);
   s.Smooth(2);
   int bm = 0;
   double mx = -1;
   for (int b = 1; b <= s.GetNbinsX(); ++b) {
      double c = s.GetBinCenter(b);
      if (c < lo || c > hi)
         continue;
      if (s.GetBinContent(b) > mx) {
         mx = s.GetBinContent(b);
         bm = b;
      }
   }
   if (!bm) {
      fwhm = -1;
      return std::nan("");
   }
   int l = bm, r = bm;
   while (l > 1 && s.GetBinContent(l) > 0.5 * mx)
      --l;
   while (r < s.GetNbinsX() && s.GetBinContent(r) > 0.5 * mx)
      ++r;
   fwhm = s.GetBinCenter(r) - s.GetBinCenter(l);
   return s.GetBinCenter(bm);
}

void dt_dvscan(TString cache = "deuteron_kin_dd.root", double dvUsed = 1.15, double dvLo = 0.80, double dvHi = 1.40,
               double dvStep = 0.05, double Ebeam = 180.0, double mEjectAmu = 2.0135532, double mResidAmu = 16.0147013,
               double exRef = 0.0, double chi2max = 5, double keMin = 3, double vzLo = 50, double vzHi = 700,
               double icLo = -1, double icHi = -1, double tcmLo = 30, double tcmHi = 55,
               TString plotOut = "plots/dt_dvscan.png")
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = 16.0147013 * u, m2 = 2.0135532 * u, m3 = mEjectAmu * u, m4 = mResidAmu * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) {
      printf("no tree pk in %s\n", cache.Data());
      return;
   }
   float ke, theta, vertexz = 0, chi2ndf, ic = -1;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &theta);
   t->SetBranchAddress("chi2ndf", &chi2ndf);
   if (t->GetBranch("vertexz"))
      t->SetBranchAddress("vertexz", &vertexz);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC)
      t->SetBranchAddress("ic", &ic);

   struct Ev {
      double pT, pz, vz;
   };
   std::vector<Ev> ev;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (chi2ndf > chi2max || ke <= 0)
         continue;
      if (hasIC && icLo > 0 && (ic < icLo || ic > icHi))
         continue;
      double p = std::sqrt((ke + m3) * (ke + m3) - m3 * m3), th = theta * TMath::DegToRad();
      ev.push_back({p * std::sin(th), p * std::cos(th), (double)vertexz});
   }
   printf("\n=== dt_dvscan: %s, %zu tracks, dv in use %.3f cm/us ===\n", gSystem->BaseName(cache), ev.size(), dvUsed);
   printf("conditions: the E_x = %.3f peak must be FLAT across theta_cm and sit AT %.3f\n", exRef, exRef);
   printf("%-8s %6s %7s %9s %8s %11s %9s\n", "dv", "k", "N", "peak", "FWHM", "slope/10deg", "spread");

   TGraph *gSlope = new TGraph(), *gMu = new TGraph(), *gW = new TGraph();
   for (double dv = dvLo; dv <= dvHi + 1e-9; dv += dvStep) {
      const double k = dv / dvUsed;
      TH1F h("h", "", 200, exRef - 6, exRef + 6);
      h.SetDirectory(nullptr);
      std::vector<std::pair<double, double>> pts; // (theta_cm, Ex)
      for (auto &e : ev) {
         double pz = k * e.pz, pT = e.pT, vz = k * e.vz;
         if (vzLo < vzHi && (vz < vzLo || vz > vzHi))
            continue;
         double p2 = pT * pT + pz * pz;
         double kep = std::sqrt(p2 + m3 * m3) - m3;
         if (kep < keMin)
            continue;
         auto [ex, tcm] = kine2bd(m1, m2, m3, m4, Ebeam, std::atan2(pT, pz), kep);
         if (std::isnan(ex))
            continue;
         h.Fill(ex);
         if (!std::isnan(tcm))
            pts.emplace_back(tcm, ex);
      }
      if (h.GetEntries() < 500) {
         printf("%-8.3f %6.3f %7.0f   (too few)\n", dv, k, h.GetEntries());
         continue;
      }
      double fw = 0, pk = peakOf(&h, fw);
      // track the peak INSIDE each theta_cm slice, fit-free, within +-1.5 MeV of the global one
      TGraph gs;
      double lo = 1e9, hi = -1e9;
      for (double a = tcmLo; a < tcmHi; a += 5) {
         TH1F hs("hs", "", 200, exRef - 6, exRef + 6);
         hs.SetDirectory(nullptr);
         for (auto &q : pts)
            if (q.first >= a && q.first < a + 5)
               hs.Fill(q.second);
         if (hs.GetEntries() < 120)
            continue;
         double fws = 0, p = peakOf(&hs, fws, pk - 1.5, pk + 1.5);
         if (std::isnan(p))
            continue;
         gs.SetPoint(gs.GetN(), a + 2.5, p);
         lo = std::min(lo, p);
         hi = std::max(hi, p);
      }
      double slope = 0;
      if (gs.GetN() >= 3) {
         TF1 l("l", "pol1", tcmLo, tcmHi);
         gs.Fit(&l, "QRN");
         slope = 10 * l.GetParameter(1);
      }
      double spread = (hi > lo) ? hi - lo : 0;
      printf("%-8.3f %6.3f %7.0f %9.3f %8.3f %11.4f %9.3f\n", dv, k, h.GetEntries(), pk, fw, slope, spread);
      gSlope->SetPoint(gSlope->GetN(), dv, slope);
      gMu->SetPoint(gMu->GetN(), dv, pk);
      gW->SetPoint(gW->GetN(), dv, fw);
   }

   auto zeroOf = [](TGraph *g, double target) -> double {
      for (int i = 0; i + 1 < g->GetN(); ++i) {
         double x0, y0, x1, y1;
         g->GetPoint(i, x0, y0);
         g->GetPoint(i + 1, x1, y1);
         if ((y0 - target) * (y1 - target) < 0)
            return x0 + (x1 - x0) * (target - y0) / (y1 - y0);
      }
      return std::nan("");
   };
   double dvFlat = zeroOf(gSlope, 0.0), dvZero = zeroOf(gMu, exRef);
   printf("\ndv making E_x flat in theta_cm : %s\n",
          std::isnan(dvFlat) ? "not bracketed in this range" : Form("%.3f cm/us", dvFlat));
   printf("dv putting the peak at %.3f     : %s\n", exRef,
          std::isnan(dvZero) ? "not bracketed in this range" : Form("%.3f cm/us", dvZero));
   if (!std::isnan(dvFlat) && !std::isnan(dvZero))
      printf("  -> the two conditions %s (difference %.3f cm/us)\n",
             std::fabs(dvFlat - dvZero) < 0.06 ? "AGREE: a drift-velocity error explains this"
                                               : "DISAGREE: dv alone does NOT set this scale",
             dvFlat - dvZero);

   auto *cv = new TCanvas("cvdv", "dt_dvscan", 1500, 480);
   cv->Divide(3, 1);
   cv->cd(1);
   gSlope->SetMarkerStyle(20);
   gSlope->SetTitle("E_{x} vs #theta_{cm} slope;drift velocity [cm/#mus];slope [MeV/10#circ]");
   gSlope->Draw("APL");
   {
      auto *z = new TLine(dvLo, 0, dvHi, 0);
      z->SetLineColor(kRed);
      z->SetLineStyle(2);
      z->Draw();
   }
   cv->cd(2);
   gMu->SetMarkerStyle(20);
   gMu->SetTitle("peak position;drift velocity [cm/#mus];E_{x} [MeV]");
   gMu->Draw("APL");
   {
      auto *z = new TLine(dvLo, exRef, dvHi, exRef);
      z->SetLineColor(kRed);
      z->SetLineStyle(2);
      z->Draw();
   }
   cv->cd(3);
   gW->SetMarkerStyle(20);
   gW->SetTitle("peak FWHM;drift velocity [cm/#mus];FWHM [MeV]");
   gW->Draw("APL");
   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
