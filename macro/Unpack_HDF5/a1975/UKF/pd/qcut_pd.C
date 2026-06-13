/// @file qcut_pd.C
/// @brief Does a FIT-QUALITY cut improve 15C(p,d) doublet resolution?
///
/// Reads deuteron_kin.root, applies the calibration + theta_cm correction (pol2
/// fit on the full in-window sample, same as excorr2_pd.C), then for a series of
/// chi2/ndf cuts measures the g.s./0.74 doublet FWHM, peak/valley, and surviving
/// yield. Prints a table so we can see resolution-vs-yield as the cut tightens.
///
///   root -b -q 'pd/qcut_pd.C'

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
static void dres(std::vector<float> &v, double &fwhm, double &ptv)
{
   if (v.size() < 200) {
      fwhm = -1;
      ptv = -1;
      return;
   }
   TH1F *h = new TH1F("ht", "", 240, -6, 14);
   h->SetDirectory(nullptr);
   for (float e : v)
      h->Fill(e);
   double gs = gsPos(h);
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

void qcut_pd(double exShift = -0.615, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ex, thcm, vz, c2n;
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("thcm", &thcm);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("chi2ndf", &c2n);
   std::vector<float> e0, tc, q;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18 && thcm > 10 && thcm < 160 && vz > 150 && vz < 1000) {
         e0.push_back(ex + exShift);
         tc.push_back(thcm);
         q.push_back(c2n);
      }
   }
   f->Close();

   // --- theta_cm correction (pol2 on g.s. position per thcm bin), fit on FULL sample ---
   TF1 sysT("sysT", "pol2", 20, 150);
   {
      TGraph *g = new TGraph();
      int nb = 13;
      double vlo = 20, vhi = 150, bw = (vhi - vlo) / nb;
      for (int b = 0; b < nb; ++b) {
         double c0 = vlo + b * bw, c1 = c0 + bw;
         TH1F *hb = new TH1F("hbb", "", 120, -3, 5);
         hb->SetDirectory(nullptr);
         long n = 0;
         for (size_t i = 0; i < e0.size(); ++i)
            if (tc[i] >= c0 && tc[i] < c1) {
               hb->Fill(e0[i]);
               ++n;
            }
         if (n > 250)
            g->SetPoint(g->GetN(), 0.5 * (c0 + c1), gsPos(hb));
         delete hb;
      }
      g->Fit(&sysT, "QRN");
   }
   // theta_cm-corrected Ex for every track
   std::vector<float> ec(e0.size());
   for (size_t i = 0; i < e0.size(); ++i) {
      double a = std::min(150.0, std::max(20.0, (double)tc[i]));
      ec[i] = e0[i] - sysT.Eval(a);
   }

   // --- scan chi2/ndf cut ---
   double cuts[] = {1e9, 5, 3, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01};
   printf("\n  chi2/ndf<     yield      doublet_FWHM   peak/valley\n");
   printf("  ------------------------------------------------------\n");
   for (double c : cuts) {
      std::vector<float> sub;
      for (size_t i = 0; i < ec.size(); ++i)
         if (q[i] < c)
            sub.push_back(ec[i]);
      double fw, pv;
      dres(sub, fw, pv);
      if (c > 1e8)
         printf("  (all)        %6zu       %6.3f        %.3f\n", sub.size(), fw, pv);
      else
         printf("  %-8g     %6zu       %6.3f        %.3f\n", c, sub.size(), fw, pv);
   }
   printf("\n  NOTE: chi2/ndf mean ~0.03 (error model loose) -> cut barely changes sample.\n");

   // --- overlay: self-calibrated Ex doublet region at 3 cuts ---
   auto mkH = [&](double cut, int col, const char *nm) {
      std::vector<float> sub;
      for (size_t i = 0; i < ec.size(); ++i)
         if (q[i] < cut)
            sub.push_back(ec[i]);
      TH1F *htmp = new TH1F("htmp", "", 240, -6, 14);
      htmp->SetDirectory(nullptr);
      for (float e : sub)
         htmp->Fill(e);
      double gs = gsPos(htmp);
      delete htmp;
      TH1F *h = new TH1F(nm, "", 140, -2, 5);
      h->SetDirectory(nullptr);
      for (float e : sub)
         h->Fill(e - gs);
      h->SetLineColor(col);
      h->SetLineWidth(2);
      return h;
   };
   TH1F *hA = mkH(1e9, kBlack, "hA");
   TH1F *hB = mkH(0.02, kBlue + 1, "hB");
   TH1F *hC = mkH(0.01, kRed, "hC");
   // scale tighter cuts up to the all-cut g.s. peak height for shape comparison
   double pkA = hA->GetMaximum();
   hB->Scale(pkA / std::max(1.0, hB->GetMaximum()));
   hC->Scale(pkA / std::max(1.0, hC->GetMaximum()));
   TCanvas *c = new TCanvas("c", "qcut", 800, 600);
   hA->SetTitle("^{15}C(p,d) g.s./0.74 doublet vs #chi^{2}/ndf cut;E_{x} - E_{x}^{gs} [MeV];counts (peak-normalized)");
   hA->Draw("hist");
   hB->Draw("hist same");
   hC->Draw("hist same");
   TLegend *lg = new TLegend(0.58, 0.7, 0.88, 0.88);
   lg->AddEntry(hA, "all (FWHM 0.61)", "l");
   lg->AddEntry(hB, "#chi^{2}/ndf<0.02 (0.53)", "l");
   lg->AddEntry(hC, "#chi^{2}/ndf<0.01 (0.41)", "l");
   lg->Draw();
   c->SaveAs("pd/plots/qcut.png");
   printf("saved pd/plots/qcut.png\n");
}
