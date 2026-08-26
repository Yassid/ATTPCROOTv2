/// @file elastic_wtrack_C14.C
/// @brief Elastic angular distribution with a WIDTH-TRACKING window: count in mu(theta) +- k*w(theta).
///
/// Why: the elastic peak broadens from sigma ~0.15 MeV forward to ~0.9 backward, so ANY fixed window
/// captures a different fraction of it at every angle. That is not a normalisation error, it is an
/// angle-dependent efficiency, and it is what makes L(theta) drift whichever fixed width is chosen:
///   +-1.0 MeV : good forward, loses a factor ~2 of the peak at 90-110
///   +-3.0 MeV : good at 90-110, over-counts beyond 125 where the peak is weak and the window wide
/// A window scaled to the MEASURED width holds the same fraction at every angle, so the efficiency
/// is constant and cancels in the luminosity. This is the refinement 09_extraction.tex recommends
/// ("a window tracking the measured peak width, mu +- 2.5 w, would be choice-free").
///
/// Pass 1 measures mu and w per bin: mode, then an iterated truncated centroid, then w from the
/// HALF-MAXIMUM crossings (insensitive to the tails, unlike an rms or a gaussian fit). Bins that
/// fail are interpolated from the neighbours that passed, then both are smoothed by a 3-point
/// moving average -- the locus is a property of the reconstruction, not of one bin's statistics.
/// Pass 2 simply COUNTS in mu +- k*w. No background subtraction: below the elastic peak there is
/// nothing to subtract.
void elastic_wtrack_C14(TString cache = "plots/proton_kin_cat5_s013.root",
                        TString accDir = "/mnt/f/a1954_C14_acc_catima/", Double_t kSig = 2.5,
                        Double_t cmMin = 15, Double_t cmMax = 150, Double_t dcm = 5.0,
                        Double_t zMin = -1e9, Double_t zMax = 1e9, Double_t chi2Cut = 5.0,
                        TString tag = "wtrack")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = (TNtuple *)fd->Get("pk");
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   TH1D *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   std::vector<double> mu(NB, 0), w(NB, 0), ctr(NB, 0);
   std::vector<bool> ok(NB, false);
   std::vector<TH1D *> hs(NB);
   float *v;
   for (int b = 0; b < NB; ++b) {
      hs[b] = new TH1D(Form("hx%d", b), "", 200, -6, 4);
      ctr[b] = cmMin + (b + 0.5) * dcm;
   }
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); v = t->GetArgs();
      if (v[5] > chi2Cut || v[0] <= 0 || v[2] < zMin || v[2] > zMax) continue;
      int b = (int)((v[3] - cmMin) / dcm);
      if (b >= 0 && b < NB) hs[b]->Fill(v[4]);
   }
   // ---- pass 1: locus and width -------------------------------------------------------------
   for (int b = 0; b < NB; ++b) {
      if (hs[b]->Integral() < 60) continue;
      TH1D *h = (TH1D *)hs[b]->Clone(Form("s%d", b)); h->Smooth(2);
      double pk = h->GetBinCenter(h->GetMaximumBin()), ymax = h->GetMaximum();
      // half-maximum crossings about the mode
      int bm = h->GetMaximumBin(), lo = bm, hi = bm;
      while (lo > 1 && h->GetBinContent(lo) > 0.5 * ymax) --lo;
      while (hi < h->GetNbinsX() && h->GetBinContent(hi) > 0.5 * ymax) ++hi;
      double fwhm = h->GetBinCenter(hi) - h->GetBinCenter(lo);
      if (fwhm <= 0) continue;
      double ww = fwhm / 2.355;
      // iterated truncated centroid inside +-1.2 w
      double m = pk;
      for (int it = 0; it < 5; ++it) {
         double s = 0, n = 0;
         for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double x = h->GetBinCenter(i);
            if (std::fabs(x - m) > 1.2 * ww) continue;
            s += x * h->GetBinContent(i); n += h->GetBinContent(i);
         }
         if (n > 0) m = s / n;
      }
      mu[b] = m; w[b] = ww; ok[b] = true;
   }
   // interpolate the failures, then smooth: these are properties of the reconstruction
   for (int b = 0; b < NB; ++b) if (!ok[b]) {
      int l = b, r = b;
      while (l >= 0 && !ok[l]) --l;
      while (r < NB && !ok[r]) ++r;
      if (l < 0 && r >= NB) continue;
      if (l < 0) { mu[b] = mu[r]; w[b] = w[r]; }
      else if (r >= NB) { mu[b] = mu[l]; w[b] = w[l]; }
      else { double f = (double)(b - l) / (r - l);
             mu[b] = mu[l] + f * (mu[r] - mu[l]); w[b] = w[l] + f * (w[r] - w[l]); }
   }
   std::vector<double> ms = mu, ws = w;
   for (int b = 1; b < NB - 1; ++b) { ms[b] = (mu[b-1]+mu[b]+mu[b+1])/3; ws[b] = (w[b-1]+w[b]+w[b+1])/3; }
   // ---- pass 2: count -----------------------------------------------------------------------
   TGraph gf; std::ifstream in((here + "/../fresco/outputs/p14C_el_161_dsdo.dat").Data());
   double a, bb; while (in >> a >> bb) if (a >= 10 && a <= 160) gf.SetPoint(gf.GetN(), a, bb);
   auto *ds = new TH1D("dsw", "", NB, cmMin, cmMax);
   printf("\n theta_cm    mu      w    window        counts   acc     dsdo      FRESCO      L\n");
   std::vector<double> Lv, Lt;
   for (int b = 0; b < NB; ++b) {
      double lo = ms[b] - kSig * ws[b], hi = ms[b] + kSig * ws[b];
      double y = hs[b]->Integral(hs[b]->FindBin(lo), hs[b]->FindBin(hi));
      double c = ctr[b];
      double dOm = 2 * TMath::Pi() * (std::cos((c - dcm/2)*TMath::DegToRad()) - std::cos((c + dcm/2)*TMath::DegToRad()));
      double A = acc ? acc->GetBinContent(acc->FindBin(c)) : 1;
      if (y <= 0 || A <= 0.05 || c < 18 || c > 148) continue;
      double d = y / A / dOm, L = d / gf.Eval(c);
      ds->SetBinContent(b + 1, d); ds->SetBinError(b + 1, std::sqrt(y) / A / dOm);
      printf("  %6.1f %+6.2f %6.3f  %+5.2f..%+5.2f %8.0f  %.3f %9.4g %9.4g %8.1f\n",
             c, ms[b], ws[b], lo, hi, y, A, d, gf.Eval(c), L);
      Lv.push_back(L); Lt.push_back(c);
   }
   // flatness over the region where the peak is measurable
   std::vector<double> Lsel;
   for (size_t i = 0; i < Lv.size(); ++i) if (Lt[i] >= 70 && Lt[i] <= 140) Lsel.push_back(Lv[i]);
   std::sort(Lsel.begin(), Lsel.end());
   if (!Lsel.empty())
      printf("\n  L over 70-140 deg: median %.1f, range %.1f - %.1f  (%+.0f%% / %+.0f%%)\n",
             Lsel[Lsel.size()/2], Lsel.front(), Lsel.back(),
             100*(Lsel.front()/Lsel[Lsel.size()/2]-1), 100*(Lsel.back()/Lsel[Lsel.size()/2]-1));
   double sd = 0, sf = 0;
   for (int b = 1; b <= ds->GetNbinsX(); ++b) { double x = ds->GetBinCenter(b);
      if (x < 80 || x > 140 || ds->GetBinContent(b) <= 0) continue; sd += ds->GetBinContent(b); sf += gf.Eval(x); }
   double sc = sd / sf;
   auto *c1 = new TCanvas("cw", "", 1000, 700); c1->SetLogy(); gPad->SetLeftMargin(0.12);
   ds->SetTitle(Form("^{14}C(p,p') elastic, window #mu #pm %.1fw;#theta_{cm} [deg];d#sigma/d#Omega [arb.]", kSig));
   ds->GetXaxis()->SetRangeUser(15, 150); ds->SetMarkerStyle(20); ds->Draw("PE");
   auto *g2 = new TGraph(); for (int i = 0; i < gf.GetN(); ++i) g2->SetPoint(i, gf.GetX()[i], sc*gf.GetY()[i]);
   g2->SetLineColor(kRed+1); g2->SetLineWidth(3); g2->Draw("L same");
   auto *lg = new TLegend(0.45,0.75,0.89,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(ds, Form("width-tracking, #pm%.1fw", kSig), "pe");
   lg->AddEntry(g2, "FRESCO KD03, scaled 80-140", "l"); lg->Draw();
   c1->SaveAs(here + "/plots/elastic_wtrack_" + tag + ".png");
}
