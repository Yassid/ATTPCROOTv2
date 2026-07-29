/// @file excorr_C15.C
/// @brief Flatten the Ex-vs-theta_cm dependence of the a2091 15C(p,p') elastic peak.
///
/// A discrete state must have Ex independent of theta_cm. Any residual trend is a kinematic
/// systematic (drift velocity / B field / energy-loss / vertex), and it smears the integrated
/// spectrum. This fits the elastic Ex(theta_cm) profile with a polynomial and subtracts it
/// event by event, then reports the g.s. width before/after and how well the excited group
/// sharpens — the same treatment used for a1975 15C(p,d).
///
/// It also prints the trend in theta_LAB, because a wrong BEAM ENERGY produces a specific
/// tilt there: distinguishing "wrong Ebeam" from "detector systematic" is the whole point.
///
///   root -b -q 'pp/excorr_C15.C("plots/proton_kin_5run_ukf.root",161,"ukf")'
static double omegaE(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

void excorr_C15(TString cache = "plots/proton_kin_5run_ukf.root", double Ebeam = 195.0, TString tag = "ukf",
                double chi2Cut = 5.0, int nPol = 2, double tcLo = 15, double tcHi = 75, double gsWin = 1.2)
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
   Float_t ke, th, vz, thcm, ex, c2;
   nt->SetBranchAddress("ke", &ke);
   nt->SetBranchAddress("theta", &th);
   nt->SetBranchAddress("thcm", &thcm);
   nt->SetBranchAddress("ex", &ex);
   nt->SetBranchAddress("chi2ndf", &c2);
   std::vector<double> vEx, vTc, vTh;
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      if (c2 > chi2Cut || !std::isfinite(ex) || ex < -5 || ex > 15)
         continue;
      vEx.push_back(ex);
      vTc.push_back(thcm);
      vTh.push_back(th);
   }
   f->Close();
   printf("%s @ Ebeam=%.1f : %zu tracks\n", tag.Data(), Ebeam, vEx.size());

   auto peakIn = [&](std::vector<double> &sel) { // gaussian peak of a sub-sample
      if (sel.size() < 40)
         return std::make_pair((double)NAN, (double)NAN);
      TH1D h("hh", "", 120, -3, 3);
      for (double e : sel)
         h.Fill(e);
      double pk = h.GetBinCenter(h.GetMaximumBin());
      TF1 g("g", "gaus", pk - 0.7, pk + 0.7);
      h.Fit(&g, "RQ0");
      return std::make_pair(g.GetParameter(1), g.GetParameter(2));
   };

   // ---------- elastic profile vs theta_cm ----------------------------------
   std::vector<double> all(vEx);
   auto [mu0, sg0] = peakIn(all);
   printf("  integrated g.s. : mu=%+.3f  sigma=%.3f  FWHM=%.3f MeV\n", mu0, sg0, 2.355 * sg0);

   const int NB = 12;
   double bw = (tcHi - tcLo) / NB;
   TGraphErrors *gp = new TGraphErrors();
   for (int b = 0; b < NB; ++b) {
      double c0 = tcLo + b * bw, c1 = c0 + bw;
      std::vector<double> sel;
      for (size_t i = 0; i < vEx.size(); ++i)
         if (vTc[i] >= c0 && vTc[i] < c1 && std::fabs(vEx[i] - mu0) < gsWin)
            sel.push_back(vEx[i]);
      auto [m, s] = peakIn(sel);
      if (!std::isfinite(m))
         continue;
      int n = gp->GetN();
      gp->SetPoint(n, 0.5 * (c0 + c1), m);
      gp->SetPointError(n, 0, s / std::sqrt((double)sel.size()));
   }
   if (gp->GetN() < nPol + 2) {
      printf("\033[1;31mtoo few theta_cm bins with an elastic peak (%d)\033[0m\n", gp->GetN());
      return;
   }
   TF1 *prof = new TF1("prof", Form("pol%d", nPol), tcLo, tcHi);
   gp->Fit(prof, "RQ0");
   double swing = 0, pmin = 1e9, pmax = -1e9;
   for (double t = tcLo; t <= tcHi; t += 0.5) {
      pmin = std::min(pmin, prof->Eval(t));
      pmax = std::max(pmax, prof->Eval(t));
   }
   swing = pmax - pmin;
   printf("  Ex(theta_cm) profile: pol%d, swing over %.0f-%.0f deg = \033[1;33m%.3f MeV\033[0m\n", nPol, tcLo, tcHi,
          swing);

   // linear trend in theta_lab too (the beam-energy signature)
   {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      long n = 0;
      for (size_t i = 0; i < vEx.size(); ++i) {
         if (vTh[i] < 55 || vTh[i] > 85 || std::fabs(vEx[i] - mu0) > gsWin)
            continue;
         sx += vTh[i];
         sy += vEx[i];
         sxx += vTh[i] * vTh[i];
         sxy += vTh[i] * vEx[i];
         ++n;
      }
      if (n > 50)
         printf("  dEx/dtheta_lab (55-85 deg, %ld trk) = %+.5f MeV/deg\n", n, (n * sxy - sx * sy) / (n * sxx - sx * sx));
   }

   // ---------- apply the correction ----------------------------------------
   std::vector<double> corr(vEx.size());
   for (size_t i = 0; i < vEx.size(); ++i) {
      double t = std::min(std::max(vTc[i], tcLo), tcHi); // clamp outside the fitted range
      corr[i] = vEx[i] - prof->Eval(t);
   }
   auto [mu1, sg1] = peakIn(corr);
   printf("  \033[1;32mafter correction: mu=%+.3f  sigma=%.3f  FWHM=%.3f MeV   (%.0f%% narrower)\033[0m\n", mu1, sg1,
          2.355 * sg1, 100 * (1 - sg1 / sg0));

   // excited group before/after
   auto bumpFit = [&](std::vector<double> &v, double lo, double hi) {
      TH1D h("hb", "", 240, -3, 15);
      for (double e : v)
         h.Fill(e);
      TF1 g("g", "gaus(0)+pol1(3)", lo, hi);
      g.SetParameters(0.3 * h.GetMaximum(), 0.5 * (lo + hi), 0.4, 5, 0);
      h.Fit(&g, "RQ0");
      return std::make_pair(g.GetParameter(1), std::fabs(g.GetParameter(2)));
   };
   auto [b0, bs0] = bumpFit(vEx, mu0 + 5.0, mu0 + 9.0);
   auto [b1, bs1] = bumpFit(corr, mu1 + 5.0, mu1 + 9.0);
   printf("  excited group : %.3f (FWHM %.2f) -> %.3f (FWHM %.2f)  [Ex above g.s.: %.3f -> %.3f]\n", b0, 2.355 * bs0,
          b1, 2.355 * bs1, b0 - mu0, b1 - mu1);

   // ---------- plots --------------------------------------------------------
   TH2F *h2a = new TH2F("h2a", TString::Format("%s: E_{x} vs #theta_{cm} BEFORE;#theta_{cm} [deg];E_{x} [MeV]", tag.Data()),
                        90, 0, 90, 160, -3, 5);
   TH2F *h2b = new TH2F("h2b", TString::Format("%s: E_{x} vs #theta_{cm} AFTER;#theta_{cm} [deg];E_{x} [MeV]", tag.Data()),
                        90, 0, 90, 160, -3, 5);
   TH1D *h1a = new TH1D("h1a", TString::Format("%s: E_{x};E_{x} [MeV];counts/50 keV", tag.Data()), 360, -3, 15);
   TH1D *h1b = new TH1D("h1b", "", 360, -3, 15);
   for (size_t i = 0; i < vEx.size(); ++i) {
      h2a->Fill(vTc[i], vEx[i]);
      h2b->Fill(vTc[i], corr[i]);
      h1a->Fill(vEx[i]);
      h1b->Fill(corr[i]);
   }
   TCanvas *c = new TCanvas("c", "excorr", 1500, 900);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLogz();
   h2a->Draw("colz");
   gp->SetMarkerStyle(20);
   gp->SetMarkerColor(kRed + 1);
   gp->SetLineColor(kRed + 1);
   gp->Draw("P same");
   prof->SetLineColor(kRed + 1);
   prof->Draw("same");
   c->cd(2);
   gPad->SetLogz();
   h2b->Draw("colz");
   c->cd(3);
   h1a->SetLineColor(kGray + 2);
   h1a->Draw();
   h1b->SetLineColor(kRed + 1);
   h1b->SetLineWidth(2);
   h1b->Draw("same");
   auto *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
   leg->AddEntry(h1a, Form("before (FWHM %.3f)", 2.355 * sg0), "l");
   leg->AddEntry(h1b, Form("after  (FWHM %.3f)", 2.355 * sg1), "l");
   leg->Draw();
   c->cd(4);
   gPad->SetLogy();
   h1b->GetXaxis()->SetRangeUser(-3, 15);
   h1b->DrawCopy();
   TString png = here + "/plots/excorr_C15_" + tag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
