/// @file spectrum_pp.C
/// @brief The a1975 16C(p,p') excitation-energy spectrum, from the full 84-run cache.
///
/// Four panels: the kinematic plane the spectrum is built from, the spectrum itself, the same
/// split by vertex position, and the same split by lab angle. The splits are there because two
/// things are known to vary across the data set and both are easier to see than to argue about.
///
/// THE ZERO-POINT IS NOT UNDERSTOOD. The ground state sits at about +0.19 MeV and no beam energy
/// between 186 and 202 MeV brings it to zero: over that whole range it moves only from +0.159 to
/// +0.234, while the 1.766 MeV level spacing is matched at about 195.3. One parameter cannot
/// satisfy both conditions. Beam energy loss along the vertex was tried as the explanation and
/// does NOT work: the loss produces a SLOPE in Ex against vertex z, and the measured slope in
/// these data is +0.06 MeV across the whole chamber, far too small and the wrong shape to be a
/// constant +0.19 offset.
///
/// Ex is computed at a FIXED beam energy, as the analysis does. No drift correction is applied
/// here; corrections belong downstream of a spectrum one can look at plainly.
///
///   root -b -q 'spectrum_pp.C("/path/pp_kin.root")'
///   root -b -q 'spectrum_pp.C("/path/pp_kin.root",195.3)'   // the spacing-calibrated energy

void spectrum_pp(TString cache, Double_t Ebeam = 192.0, Double_t chi2Max = 5.0, Double_t icLo = 1000,
                 Double_t icHi = 1350, Double_t thLo = 25, Double_t thHi = 65, TString tag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t)
      return;
   float ke, th, vz, c2, ic;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);
   t->SetBranchAddress("ic", &ic);

   const double u = 931.49401;
   const double mb = 16.0147 * u, mt = 1.0078250322 * u, m3 = 1.0078250322 * u, m4 = 16.0147 * u;
   auto om = [](double x, double y, double z) { return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z); };
   auto exOf = [&](double K, double thr, double Ke) {
      double E1 = K + mb, E3 = Ke + m3;
      double s = mb * mb + mt * mt + 2 * mt * E1, uu = mt * mt + m3 * m3 - 2 * mt * E3;
      double a = (std::cos(thr) * om(s, mb * mb, mt * mt) * om(uu, mt * mt, m3 * m3) -
                  (s - mb * mb - mt * mt) * (mt * mt + m3 * m3 - uu)) /
                    (2 * mt * mt) +
                 s + uu - mt * mt;
      return a < 0 ? std::nan("") : std::sqrt(a) - m4;
   };

   auto *hkt = new TH2D("hkt", "p + ^{16}C kinematics;#theta_{lab} [deg];KE_{p} [MeV]", 100, 0, 90, 100, 0, 50);
   auto *hEx = new TH1D("hEx", "", 300, -3, 9);
   auto *hz1 = new TH1D("hz1", "", 300, -3, 9);
   auto *hz2 = new TH1D("hz2", "", 300, -3, 9);
   auto *ht1 = new TH1D("ht1", "", 300, -3, 9);
   auto *ht2 = new TH1D("ht2", "", 300, -3, 9);
   long n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2Max && ic > icLo && ic < icHi))
         continue;
      hkt->Fill(th, ke);
      if (th < thLo || th > thHi)
         continue;
      double e = exOf(Ebeam, th * TMath::DegToRad(), ke);
      if (std::isnan(e))
         continue;
      hEx->Fill(e);
      ++n;
      (vz < 450 ? hz1 : hz2)->Fill(e);
      (th < 45 ? ht1 : ht2)->Fill(e);
   }
   printf("\n  %ld protons  (chi2/ndf < %.0f, IC %.0f-%.0f, theta_lab %.0f-%.0f)\n", n, chi2Max, icLo, icHi, thLo, thHi);

   // the known 16C levels, for reference lines -- NOT fitted
   const int NL = 3;
   const double LV[NL] = {0.0, 1.766, 3.027};
   const char *LN[NL] = {"g.s. 0^{+}", "1.766 (2^{+})", "3.027 (0^{+})"};

   auto peak = [&](TH1D *h, double c, double w, double &mu, double &sg) {
      TF1 g("g", "gaus(0)+pol1(3)", c - w, c + w);
      int bm = h->GetMaximumBin();
      g.SetParameters(h->GetBinContent(bm), c, 0.3, 0, 0);
      g.SetParLimits(1, c - w, c + w);
      g.SetParLimits(2, 0.05, 1.5);
      if (h->Fit(&g, "QRN") != 0) { mu = sg = -99; return false; }
      mu = g.GetParameter(1); sg = g.GetParameter(2);
      return true;
   };
   double mu0, sg0, mu1, sg1;
   bool ok0 = peak(hEx, 0.2, 0.7, mu0, sg0), ok1 = peak(hEx, 1.9, 0.7, mu1, sg1);
   if (ok0 && ok1)
      printf("  g.s. at %+.3f (sigma %.3f),  2+ at %+.3f (sigma %.3f),  spacing %.3f vs 1.766 known\n", mu0, sg0, mu1,
             sg1, mu1 - mu0);

   TCanvas *c = new TCanvas("cpp", "pp spectrum", 1500, 950);
   c->Divide(2, 2);

   c->cd(1);
   gPad->SetLogz();
   hkt->Draw("colz");

   auto marks = [&](double ymax) {
      for (int i = 0; i < NL; ++i) {
         auto *l = new TLine(LV[i], 0, LV[i], ymax);
         l->SetLineColor(kRed + 1);
         l->SetLineStyle(2);
         l->SetLineWidth(2);
         l->Draw();
         auto *tx = new TLatex(LV[i] + 0.06, ymax * 0.88, LN[i]);
         tx->SetTextColor(kRed + 1);
         tx->SetTextSize(0.033);
         tx->SetTextAngle(90);
         tx->Draw();
      }
   };

   c->cd(2);
   hEx->SetTitle(TString::Format("E_{x}(^{16}C) at E_{beam} = %.1f MeV;E_{x} [MeV];protons", Ebeam));
   hEx->SetLineColor(kBlack);
   hEx->SetLineWidth(2);
   hEx->Draw("hist");
   marks(hEx->GetMaximum() * 1.05);

   c->cd(3);
   double m3v = std::max(hz1->GetMaximum(), hz2->GetMaximum());
   hz1->SetTitle("split by vertex z;E_{x} [MeV];protons");
   hz1->SetMaximum(m3v * 1.25);
   hz1->SetLineColor(kAzure + 2);
   hz1->SetLineWidth(2);
   hz1->Draw("hist");
   hz2->SetLineColor(kRed + 1);
   hz2->SetLineWidth(2);
   hz2->Draw("hist same");
   marks(m3v * 1.25);
   auto *l3 = new TLegend(0.60, 0.72, 0.89, 0.88);
   l3->AddEntry(hz1, "z < 450 mm", "l");
   l3->AddEntry(hz2, "z > 450 mm", "l");
   l3->Draw();

   c->cd(4);
   double m4v = std::max(ht1->GetMaximum(), ht2->GetMaximum());
   ht1->SetTitle("split by lab angle;E_{x} [MeV];protons");
   ht1->SetMaximum(m4v * 1.25);
   ht1->SetLineColor(kGreen + 3);
   ht1->SetLineWidth(2);
   ht1->Draw("hist");
   ht2->SetLineColor(kOrange + 7);
   ht2->SetLineWidth(2);
   ht2->Draw("hist same");
   marks(m4v * 1.25);
   auto *l4 = new TLegend(0.60, 0.72, 0.89, 0.88);
   l4->AddEntry(ht1, TString::Format("#theta_{lab} < 45#circ"), "l");
   l4->AddEntry(ht2, TString::Format("#theta_{lab} > 45#circ"), "l");
   l4->Draw();

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/spectrum_pp" + tag + ".png";
   c->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());
}
