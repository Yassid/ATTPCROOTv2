/// @file excorr2_pd.C
/// @brief Combined theta_cm + vertex-Z systematic correction for 15C(p,d) Ex.
///
/// Both theta_cm (U-shape) and vertex-Z (rising) introduce a g.s.-position drift
/// that smears the integrated spectrum. This flattens them sequentially: per-theta_cm
/// bin g.s. position -> pol2 subtract; then per-Z bin g.s. position on the result
/// -> pol1 subtract. The g.s. peak position each bin is found with the STABLE
/// fixed-0.74-sep common-sigma doublet fit. Reports the doublet resolution gain.
///
///   root -b -q 'pd/excorr2_pd.C'

// g.s. position via fixed-0.74-sep, common-sigma, linear-bg doublet fit (stable)
static double gsPos(TH1 *h)
{
   double mx = h->GetMaximum();
   TF1 dg("dg", "[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x", -1.5, 2.3);
   dg.SetParameters(0.45, 0.32, 0.6 * mx, 0.6 * mx, 0.1 * mx, 0.0);
   dg.SetParLimits(0, -0.6, 1.2);
   dg.SetParLimits(1, 0.12, 0.55);
   dg.SetParLimits(2, 0, 2 * mx);
   dg.SetParLimits(3, 0, 2 * mx);
   h->Fit(&dg, "QRN");
   return dg.GetParameter(0);
}
// doublet resolution (FWHM, PtV) of a self-calibrated spectrum
static void dres(std::vector<float> &v, double &fwhm, double &ptv, double &gsout)
{
   TH1F *h = new TH1F("ht", "", 240, -6, 14);
   h->SetDirectory(nullptr);
   for (float e : v)
      h->Fill(e);
   double gs = gsPos(h);
   gsout = gs;
   double mx = h->GetMaximum();
   TF1 dg("dg2", "[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x", gs - 1.5, gs + 1.6);
   dg.SetParameters(gs, 0.30, 0.6 * mx, 0.6 * mx, 0.1 * mx, 0.0);
   dg.SetParLimits(0, gs - 0.4, gs + 0.4);
   dg.SetParLimits(1, 0.12, 0.55);
   h->Fit(&dg, "QRN");
   fwhm = 2.3548 * dg.GetParameter(1);
   TH1F *hs = new TH1F("hts", "", 240, -6 - gs, 14 - gs);
   hs->SetDirectory(nullptr);
   for (float e : v)
      hs->Fill(e - gs);
   hs->Smooth(1);
   auto mm = [&](double lo, double hi, bool mn) {
      double m = mn ? 1e18 : 0;
      for (int b = hs->FindBin(lo); b <= hs->FindBin(hi); ++b)
         m = mn ? std::min(m, hs->GetBinContent(b)) : std::max(m, hs->GetBinContent(b));
      return m;
   };
   ptv = std::min(mm(-0.25, 0.25, false), mm(0.5, 1.0, false)) / std::max(1.0, mm(0.25, 0.55, true));
   delete h;
   delete hs;
}

void excorr2_pd(double exShift = -0.615, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ex, thcm, vz;
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("thcm", &thcm);
   t->SetBranchAddress("vertexz", &vz);
   std::vector<float> e0, tc, zz;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18 && thcm > 10 && thcm < 160 && vz > 150 && vz < 1000) {
         e0.push_back(ex + exShift);
         tc.push_back(thcm);
         zz.push_back(vz);
      }
   }
   f->Close();
   printf("tracks: %zu\n", e0.size());

   // --- per-theta_cm bin g.s. position -> pol2 ---
   auto fitVar = [&](std::vector<float> &E, std::vector<float> &V, double vlo, double vhi, int nb, const char *form,
                     TF1 &out) {
      TGraph *g = new TGraph();
      double bw = (vhi - vlo) / nb;
      for (int b = 0; b < nb; ++b) {
         double c0 = vlo + b * bw, c1 = c0 + bw;
         TH1F *hb = new TH1F("hbb", "", 120, -3, 5);
         hb->SetDirectory(nullptr);
         long n = 0;
         for (size_t i = 0; i < E.size(); ++i)
            if (V[i] >= c0 && V[i] < c1) {
               hb->Fill(E[i]);
               ++n;
            }
         if (n > 250)
            g->SetPoint(g->GetN(), 0.5 * (c0 + c1), gsPos(hb));
         delete hb;
      }
      g->Fit(&out, "QRN");
      return g;
   };

   TF1 sysT("sysT", "pol2", 20, 150);
   TGraph *gT = fitVar(e0, tc, 20, 150, 13, "pol2", sysT);
   // apply theta_cm correction
   std::vector<float> e1(e0.size());
   for (size_t i = 0; i < e0.size(); ++i) {
      double a = std::min(150.0, std::max(20.0, (double)tc[i]));
      e1[i] = e0[i] - sysT.Eval(a);
   }
   TF1 sysZ("sysZ", "pol1", 200, 1000);
   TGraph *gZ = fitVar(e1, zz, 200, 1000, 12, "pol1", sysZ);
   // apply Z correction
   std::vector<float> e2(e1.size());
   for (size_t i = 0; i < e1.size(); ++i) {
      double z = std::min(1000.0, std::max(200.0, (double)zz[i]));
      e2[i] = e1[i] - sysZ.Eval(z);
   }

   double f0, p0, g0, f1, p1, g1, f2, p2, g2;
   dres(e0, f0, p0, g0);
   dres(e1, f1, p1, g1);
   dres(e2, f2, p2, g2);
   printf("\n                       doublet_FWHM   peak/valley\n");
   printf("  raw (calibrated)        %6.3f         %.3f\n", f0, p0);
   printf("  + theta_cm corrected    %6.3f         %.3f\n", f1, p1);
   printf("  + theta_cm + Z corr     %6.3f         %.3f\n", f2, p2);
   printf("  => FWHM %.0f%% better, peak/valley %+.0f%% \n", 100 * (f0 - f2) / f0, 100 * (p2 - p0) / p0);

   // draw: corrections + before/after spectra
   TH1F *hb0 = new TH1F("hb0", "", 200, -3, 8);
   hb0->SetDirectory(nullptr);
   TH1F *hb2 = new TH1F("hb2", "", 200, -3, 8);
   hb2->SetDirectory(nullptr);
   for (size_t i = 0; i < e0.size(); ++i) {
      hb0->Fill(e0[i] - g0);
      hb2->Fill(e2[i] - g2);
   }
   TCanvas *c = new TCanvas("c", "excorr2", 1500, 520);
   c->Divide(3, 1);
   c->cd(1);
   gT->SetTitle("g.s. vs #theta_{cm};#theta_{cm};E_{x}^{gs}");
   gT->SetMarkerStyle(20);
   gT->Draw("AP");
   sysT.SetLineColor(kRed);
   sysT.Draw("same");
   c->cd(2);
   gZ->SetTitle("g.s. vs Z (after #theta_{cm});Z [mm];E_{x}^{gs}");
   gZ->SetMarkerStyle(20);
   gZ->Draw("AP");
   sysZ.SetLineColor(kRed);
   sysZ.Draw("same");
   c->cd(3);
   hb0->SetLineColor(kGray + 2);
   hb0->SetLineWidth(2);
   hb0->SetTitle("Ex: raw (grey) vs corrected (blue);E_{x}(^{15}C) [MeV];counts");
   hb0->Draw("hist");
   hb2->SetLineColor(kBlue + 1);
   hb2->SetLineWidth(2);
   hb2->Draw("hist same");
   c->SaveAs("pd/plots/excorr2.png");
   printf("saved pd/plots/excorr2.png\n");
}
