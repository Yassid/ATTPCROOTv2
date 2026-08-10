/// @file gs_cut_pp.C
/// @brief The Ex vs theta_cm plane with the tuned correction, to place a ground-state cut on it.
///
/// Applies the parameters saved from the browser explorer, so the plane here is the one that was
/// tuned there rather than an approximation of it:
///
///     Ebeam    = 188 MeV          (not the nominal 192)
///     theta    -> theta - slope*(KE - pivot),  slope = 360/kcDenom deg/MeV, pivot = kcPivot
///     kcDenom  = 4000  ->  slope = 0.090 deg/MeV,  pivot = 1.5 MeV
///     chi2/ndf < 5
///
/// The correction is applied to theta BEFORE Ex and theta_cm are computed, which is what the
/// explorer does: thCorr() wraps every use of the angle, not just the display panel.
///
/// The ground state does not sit at a constant Ex across the angular range, so a horizontal band
/// is the wrong cut: it would take a different part of the peak at each angle and hand the
/// angular distribution a selection efficiency that varies with angle -- exactly the kind of
/// artefact that looks like physics. This measures the g.s. ridge bin by bin in theta_cm and puts
/// the cut a fixed number of widths around THAT, so the same fraction of the peak is taken
/// everywhere.
///
///   root -b -q 'gs_cut_pp.C("/path/pp_kin.root")'
///   root -b -q 'gs_cut_pp.C("/path/pp_kin.root",188,4000,1.5,2.0)'   // 2 sigma instead of 2.5

void gs_cut_pp(TString cache, Double_t Ebeam = 188.0, Double_t kcDenom = 4000.0, Double_t kcPivot = 1.5,
               Double_t nSig = 2.5, Double_t chi2Max = 5.0, Double_t icLo = 1000, Double_t icHi = 1350,
               Double_t cmLo = 0, Double_t cmHi = 180, Double_t dcm = 5.0, Double_t searchWin = 0.7,
               TString tag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double slope = 360.0 / kcDenom; // deg/MeV, exactly as kcSlopeDeg() in the explorer

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
   const double mb = 16.0147013 * u, mt = 1.00782503 * u, m3 = 1.00782503 * u, m4 = 16.0147013 * u;
   auto om = [](double x, double y, double z) {
      return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
   };
   // Ex and theta_cm together, from the (corrected) lab angle
   auto kine = [&](double K, double thr, double Ke) {
      double E1 = K + mb, E3 = Ke + m3, E4 = E1 + mt - E3;
      double s = mb * mb + mt * mt + 2 * mt * E1, uu = mt * mt + m3 * m3 - 2 * mt * E3;
      double a = (std::cos(thr) * om(s, mb * mb, mt * mt) * om(uu, mt * mt, m3 * m3) -
                  (s - mb * mb - mt * mt) * (mt * mt + m3 * m3 - uu)) /
                    (2 * mt * mt) +
                 s + uu - mt * mt;
      if (a < 0)
         return std::make_pair(std::nan(""), std::nan(""));
      double m4x = std::sqrt(a), ex = m4x - m4;
      double tt = mt * mt + m4x * m4x - 2 * mt * E4;
      double arg = (s * s + s * (2 * tt - mb * mb - mt * mt - m3 * m3 - m4x * m4x) +
                    (mb * mb - mt * mt) * (m3 * m3 - m4x * m4x)) /
                   (om(s, mb * mb, mt * mt) * om(s, m3 * m3, m4x * m4x));
      if (arg < -1 || arg > 1)
         return std::make_pair(ex, std::nan(""));
      return std::make_pair(ex, (TMath::Pi() - std::acos(arg)) * TMath::RadToDeg());
   };

   const int NB = (int)std::lround((cmHi - cmLo) / dcm);
   auto *h2 = new TH2D("h2", "E_{x} vs #theta_{cm} (corrected);#theta_{cm} [deg];E_{x} [MeV]", NB, cmLo, cmHi, 300,
                       -3, 6);
   std::vector<TH1D *> sl(NB);
   for (int b = 0; b < NB; ++b)
      sl[b] = new TH1D(TString::Format("sl%d", b), "", 200, -2.5, 2.5);

   long n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2Max && ic > icLo && ic < icHi))
         continue;
      double thc = th - slope * (ke - kcPivot); // the explorer's thCorrForce
      auto [ex, cm] = kine(Ebeam, thc * TMath::DegToRad(), ke);
      if (std::isnan(ex) || std::isnan(cm))
         continue;
      h2->Fill(cm, ex);
      int b = (int)((cm - cmLo) / dcm);
      if (b >= 0 && b < NB && std::fabs(ex) < 2.5)
         sl[b]->Fill(ex);
      ++n;
   }
   printf("\n  %ld protons after chi2 and IC cuts\n", n);
   printf("  correction: Ebeam %.1f, theta -> theta - %.4f*(KE - %.2f)\n\n", Ebeam, slope, kcPivot);

   // the g.s. ridge, bin by bin
   auto *gMu = new TGraph();
   auto *gLo = new TGraph();
   auto *gHi = new TGraph();
   int m = 0;
   printf("  theta_cm |    n   |  mu     sigma  |  cut window\n");
   for (int b = 0; b < NB; ++b) {
      double c = cmLo + (b + 0.5) * dcm;
      if (sl[b]->Integral() < 150)
         continue;
      // Seed from the tallest bin WITHIN searchWin of zero, not the tallest bin overall. Where
      // the ground state is weak -- around theta_cm 50-65 and beyond 130 -- the 1.766 MeV 2+ is
      // the tallest thing in the slice, and an unconstrained seed locks onto it and reports the
      // 2+ ridge as the g.s. one.
      int bm = -1;
      double best = -1;
      for (int k = 1; k <= sl[b]->GetNbinsX(); ++k) {
         if (std::fabs(sl[b]->GetBinCenter(k)) > searchWin)
            continue;
         if (sl[b]->GetBinContent(k) > best) {
            best = sl[b]->GetBinContent(k);
            bm = k;
         }
      }
      if (bm < 0)
         continue;
      double pk = sl[b]->GetBinCenter(bm);
      TF1 g("g", "gaus(0)+pol1(3)", pk - 0.6, pk + 0.6);
      g.SetParameters(sl[b]->GetBinContent(bm), pk, 0.25, 0, 0);
      g.SetParLimits(1, -searchWin, searchWin); // never let it wander onto the 2+
      g.SetParLimits(2, 0.05, 1.2);
      if (sl[b]->Fit(&g, "QRN") != 0)
         continue;
      double mu = g.GetParameter(1), sg = g.GetParameter(2), er = g.GetParError(1);
      if (!(er > 0 && er < 0.2))
         continue; // a fit that did not converge is not a ridge point
      gMu->SetPoint(m, c, mu);
      gLo->SetPoint(m, c, mu - nSig * sg);
      gHi->SetPoint(m, c, mu + nSig * sg);
      ++m;
      printf("   %3.0f-%3.0f | %6.0f | %+6.3f  %5.3f | %+6.3f .. %+6.3f\n", c - dcm / 2, c + dcm / 2,
             sl[b]->Integral(), mu, sg, mu - nSig * sg, mu + nSig * sg);
   }
   printf("\n  %d usable angular bins\n", m);

   TCanvas *c1 = new TCanvas("cgs", "g.s. selection", 1400, 600);
   c1->Divide(2, 1);
   c1->cd(1);
   gPad->SetLogz();
   h2->Draw("colz");
   auto style = [](TGraph *g, int col, int w, int st) {
      g->SetLineColor(col);
      g->SetLineWidth(w);
      g->SetLineStyle(st);
      g->Draw("L same");
   };
   style(gMu, kRed + 1, 3, 1);
   style(gLo, kRed + 1, 3, 2);
   style(gHi, kRed + 1, 3, 2);

   c1->cd(2);
   auto *hz = (TH2D *)h2->Clone("hz");
   hz->SetTitle("zoom on the ground state");
   hz->GetYaxis()->SetRangeUser(-1.5, 1.5);
   gPad->SetLogz();
   hz->Draw("colz");
   style(gMu, kRed + 1, 3, 1);
   style(gLo, kRed + 1, 3, 2);
   style(gHi, kRed + 1, 3, 2);

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/gs_cut_pp" + tag + ".png";
   c1->SaveAs(png);

   // store the band so the angular-distribution step uses exactly this selection
   TFile fo(here + "/plots/gs_cut_pp" + tag + ".root", "RECREATE");
   gMu->Write("gs_mu");
   gLo->Write("gs_lo");
   gHi->Write("gs_hi");
   h2->Write("ex_vs_thcm");
   TNamed par(TString("params"),
              TString::Format("Ebeam=%.1f slope=%.4f pivot=%.2f nSig=%.1f chi2=%.1f", Ebeam, slope, kcPivot, nSig,
                              chi2Max));
   par.Write();
   fo.Close();
   printf("\n  wrote %s\n         plots/gs_cut_pp%s.root (band + params)\n\n", png.Data(), tag.Data());
}
