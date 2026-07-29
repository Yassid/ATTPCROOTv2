/// @file kin_C15.C
/// @brief Kinematics check for a2091 15C(p,p'): the measured recoil-proton KE-vs-theta_lab
///        locus with the two-body loci of the known 15C levels drawn on top, the Ex spectrum,
///        and a Newton solve for the beam energy that puts the elastic peak at Ex = 0.
///
///        Runs on the ntuple cache written by pp/ex_C15.C (ke:theta:vertexz:thcm:ex:chi2ndf),
///        so it costs nothing to re-run while more runs finish.
///
///   root -b -q 'pp/kin_C15.C("plots/proton_kin_first_ukf.root",161,"ukf")'

static TGraph *kinLine(double Eb, double Ex, double m1, double m2, double m3, double m4_0, Color_t col, int style)
{
   double m4 = m4_0 + Ex;
   double E1 = Eb + m1, P = std::sqrt(E1 * E1 - m1 * m1), W = E1 + m2, s = W * W - P * P;
   auto *g = new TGraph();
   g->SetLineColor(col);
   g->SetLineWidth(2);
   g->SetLineStyle(style);
   if (s <= (m3 + m4) * (m3 + m4))
      return g; // below threshold
   double rs = std::sqrt(s), E3cm = (s + m3 * m3 - m4 * m4) / (2 * rs);
   double p3cm = std::sqrt(std::max(0., E3cm * E3cm - m3 * m3));
   double beta = P / W, gamma = W / rs;
   for (double tc = 0; tc <= 180; tc += 0.25) {
      double c = std::cos(tc * TMath::DegToRad()), sn = std::sin(tc * TMath::DegToRad());
      double Elab = gamma * (E3cm + beta * p3cm * c), pz = gamma * (p3cm * c + beta * E3cm), pperp = p3cm * sn;
      double th = std::atan2(pperp, pz) * TMath::RadToDeg(), ke = Elab - m3;
      if (ke > 0 && th >= 0)
         g->SetPoint(g->GetN(), th, ke);
   }
   return g;
}

static double omega_(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

// Ex from (theta_lab, KE) for m1(m2,m3)m4 -- same expression as ex_C15.C
static double exOf(double m1, double m2, double m3, double m4, double Kp, double thRad, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thRad) * omega_(s, m1 * m1, m2 * m2) * omega_(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   return arg > 0 ? std::sqrt(arg) - m4 : NAN;
}

void kin_C15(TString cache = "plots/proton_kin_first_ukf.root", double Ebeam = 195.0, TString tag = "ukf",
             double chi2Cut = 5.0, double mEjectAmu = 1.007825, double mResidAmu = 15.0105993)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(__FILE__);
   if (!cache.BeginsWith("/"))
      cache = here + "/" + cache;
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mERROR: cannot open %s\033[0m\n", cache.Data());
      return;
   }
   TNtuple *nt = (TNtuple *)f->Get("pk");
   const double u = 931.49401;
   const double m1 = 15.0105993 * u, m2 = 1.007825 * u, m3 = mEjectAmu * u, m4 = mResidAmu * u;

   // pull the sample once (Ex is recomputed here so the beam energy can be scanned)
   Float_t ke, th, vz, thcm, ex, c2;
   nt->SetBranchAddress("ke", &ke);
   nt->SetBranchAddress("theta", &th);
   nt->SetBranchAddress("chi2ndf", &c2);
   std::vector<double> vKe, vTh;
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0 || th < 0)
         continue;
      vKe.push_back(ke);
      vTh.push_back(th);
   }
   printf("%s: %zu tracks (chi2/ndf<%g)\n", tag.Data(), vKe.size(), chi2Cut);

   // g.s. centroid as a function of Ebeam -> Newton solve for Ex(g.s.) = 0
   auto gsMean = [&](double Eb) {
      TH1D h("hgs", "", 300, -3, 3);
      for (size_t i = 0; i < vKe.size(); ++i) {
         double e = exOf(m1, m2, m3, m4, Eb, vTh[i] * TMath::DegToRad(), vKe[i]);
         if (std::isfinite(e))
            h.Fill(e);
      }
      int bmax = h.GetMaximumBin();
      TF1 g("g", "gaus", h.GetBinCenter(bmax) - 0.6, h.GetBinCenter(bmax) + 0.6);
      h.Fit(&g, "RQ0");
      return std::make_pair(g.GetParameter(1), g.GetParameter(2));
   };
   auto [mu0, sg0] = gsMean(Ebeam);
   printf("  at Ebeam=%.2f MeV : g.s. mu = %+.3f MeV, sigma = %.3f (FWHM %.3f)\n", Ebeam, mu0, sg0, 2.355 * sg0);
   double e1 = Ebeam, m_1 = mu0, e2 = Ebeam + 4.0, m_2 = gsMean(e2).first, eSol = Ebeam;
   for (int it = 0; it < 15; ++it) {
      if (std::fabs(m_2 - m_1) < 1e-9)
         break;
      double e3 = e2 - m_2 * (e2 - e1) / (m_2 - m_1);
      if (!std::isfinite(e3) || e3 < 20 || e3 > 600)
         break;
      e1 = e2;
      m_1 = m_2;
      e2 = e3;
      m_2 = gsMean(e2).first;
      if (std::fabs(m_2) < 1e-4)
         break;
   }
   eSol = e2;
   auto [muS, sgS] = gsMean(eSol);
   printf("  \033[1;32mzero-g.s. solution: Ebeam = %.2f MeV (%.3f MeV/u) -> mu %+.4f, sigma %.3f (FWHM %.3f)\033[0m\n",
          eSol, eSol / 14.0, muS, sgS, 2.355 * sgS);
   printf("  dEx/dEbeam = %.4f MeV/MeV  <- the centroid alone constrains Ebeam only to ~+-%.0f MeV\n",
          (gsMean(Ebeam + 5).first - mu0) / 5.0, std::fabs(0.05 / std::max(1e-6, (gsMean(Ebeam + 5).first - mu0) / 5.0)));

   // ---- the estimator that DOES constrain Ebeam: Ex of the elastic ridge must not depend
   //      on the scattering angle. A wrong beam energy tilts Ex vs theta_lab; solve slope = 0.
   auto exSlope = [&](double Eb) {
      auto [mu, sg] = gsMean(Eb);
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      long n = 0;
      for (size_t i = 0; i < vKe.size(); ++i) {
         if (vTh[i] < 55 || vTh[i] > 85) // the clean elastic arc
            continue;
         double e = exOf(m1, m2, m3, m4, Eb, vTh[i] * TMath::DegToRad(), vKe[i]);
         if (!std::isfinite(e) || std::fabs(e - mu) > 3 * sg)
            continue;
         sx += vTh[i];
         sy += e;
         sxx += vTh[i] * vTh[i];
         sxy += vTh[i] * e;
         ++n;
      }
      if (n < 50)
         return std::make_pair((double)NAN, (long)n);
      double den = n * sxx - sx * sx;
      return std::make_pair((n * sxy - sx * sy) / den, n); // dEx/dtheta [MeV/deg]
   };
   auto [sl0, nSl] = exSlope(Ebeam);
   printf("  dEx/dtheta_lab at Ebeam=%.1f : %+.5f MeV/deg  (%ld elastic tracks, 55-85 deg)\n", Ebeam, sl0, nSl);
   double b1 = Ebeam, s1 = sl0, b2 = Ebeam * 1.15, s2 = exSlope(b2).first, bSol = NAN;
   for (int it = 0; it < 20 && std::isfinite(s1) && std::isfinite(s2); ++it) {
      if (std::fabs(s2 - s1) < 1e-12)
         break;
      double b3 = b2 - s2 * (b2 - b1) / (s2 - s1);
      if (!std::isfinite(b3) || b3 < 20 || b3 > 800)
         break;
      b1 = b2;
      s1 = s2;
      b2 = b3;
      s2 = exSlope(b2).first;
      if (std::fabs(s2) < 1e-6)
         break;
   }
   if (std::isfinite(s2) && std::fabs(s2) < 1e-4) {
      bSol = b2;
      auto [muF, sgF] = gsMean(bSol);
      printf("  \033[1;36mflat-Ex(theta) solution: Ebeam = %.2f MeV (%.3f MeV/u) -> slope %+.6f, g.s. mu %+.3f, "
             "sigma %.3f (FWHM %.3f)\033[0m\n",
             bSol, bSol / 14.0, s2, muF, sgF, 2.355 * sgF);
   } else {
      printf("  flat-Ex(theta) solve did not converge (last slope %+.6f at %.1f MeV)\n", s2, b2);
   }

   // ---- plots at the INPUT Ebeam --------------------------------------------
   TH2F *hkt = new TH2F("hkt", TString::Format("a2091 15C(p,p') %s: recoil-proton KE vs #theta_{lab} (E_{beam}=%.1f "
                                               "MeV);#theta_{lab} [deg];KE [MeV]",
                                               tag.Data(), Ebeam),
                        180, 0, 90, 200, 0, 45);
   TH1F *hex = new TH1F("hex", TString::Format("%s: ^{14}C excitation energy;E_{x} [MeV];counts/50 keV", tag.Data()),
                        400, -5, 15);
   for (size_t i = 0; i < vKe.size(); ++i) {
      hkt->Fill(vTh[i], vKe[i]);
      double e = exOf(m1, m2, m3, m4, Ebeam, vTh[i] * TMath::DegToRad(), vKe[i]);
      if (std::isfinite(e))
         hex->Fill(e);
   }

   TCanvas *c = new TCanvas("c", "kin", 1500, 620);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz();
   hkt->Draw("colz");
   struct L {
      double ex;
      Color_t col;
      int st;
      const char *lab;
   };
   std::vector<L> lev = {{0.0, kRed + 1, 1, "g.s."},
                         {6.09, kGreen + 2, 2, "6.09 (1^{-})"},
                         {6.59, kAzure + 1, 2, "6.59 (0^{+})"},
                         {7.34, kMagenta + 1, 2, "7.34 (2^{+})"}};
   auto *leg = new TLegend(0.55, 0.62, 0.88, 0.88);
   leg->SetHeader("two-body loci");
   for (auto &l : lev) {
      TGraph *g = kinLine(Ebeam, l.ex, m1, m2, m3, m4, l.col, l.st);
      if (g->GetN()) {
         g->Draw("L same");
         leg->AddEntry(g, l.lab, "l");
      }
   }
   leg->Draw();
   c->cd(2);
   hex->Draw();
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.035);
   tx->DrawLatex(0.5, 0.85, Form("g.s. #mu=%+.3f #sigma=%.3f", mu0, sg0));
   tx->DrawLatex(0.5, 0.80, Form("zero-g.s. E_{beam}=%.1f MeV", eSol));
   TString png = here + "/plots/kin_C15_" + tag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
