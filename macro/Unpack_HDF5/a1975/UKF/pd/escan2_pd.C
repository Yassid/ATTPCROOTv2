/// @file escan2_pd.C
/// @brief Beam-energy scan (100 keV steps) for best 15C(p,d) g.s./0.74 doublet
///        resolution — STABLE common-sigma doublet-fit metric.
///
/// For each Ebeam, recomputes Ex from cached (ke, theta), fits the g.s./0.74 doublet
/// with separation FIXED 0.74, COMMON sigma, linear bg. The doublet FWHM (= 2.355*sigma)
/// integrated over all angles minimizes at the Ebeam that flattens the theta_cm
/// dependence. Reports FWHM and g.s. position vs Ebeam; picks the minimum.
///
///   root -b -q 'pd/escan2_pd.C(189.0,197.0,0.1)'

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static const double U = 931.49401, M1 = 16.0147 * U, M2 = 1.007825 * U, M3 = 2.01410178 * U, M4 = 15.0105993 * U;
static double exOf(double Kp, double th, double Ke)
{
   double Et1 = Kp + M1, Et3 = Ke + M3, s = M1 * M1 + M2 * M2 + 2 * M2 * Et1, uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double m4e = std::sqrt((std::cos(th) * omega2(s, M1 * M1, M2 * M2) * omega2(uu, M2 * M2, M3 * M3) -
                           (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                            (2 * M2 * M2) +
                          s + uu - M2 * M2);
   return m4e - M4;
}

void escan2_pd(double Elo = 189.0, double Ehi = 197.0, double dstep = 0.1, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ke, th, ex;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("ex", &ex);
   std::vector<float> vke, vth;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18 && ke > 0 && ke < 200 && th > 0 && th < 90) {
         vke.push_back(ke);
         vth.push_back(th);
      }
   }
   f->Close();

   int nE = (int)std::round((Ehi - Elo) / dstep) + 1;
   TGraph *gF = new TGraph();
   printf("  Ebeam   gs_pos   doublet_FWHM\n");
   double bestF = 1e9, bestE = Elo;
   for (int ie = 0; ie < nE; ++ie) {
      double E = Elo + ie * dstep;
      TH1F *h = new TH1F("h", "", 200, -6, 14);
      h->SetDirectory(nullptr);
      for (size_t i = 0; i < vke.size(); ++i) {
         double e = exOf(E, vth[i] * TMath::DegToRad(), vke[i]);
         if (!std::isnan(e))
            h->Fill(e);
      }
      // locate doublet, then common-sigma fixed-0.74 fit on a window around it
      h->GetXaxis()->SetRangeUser(-1.5, 2.0);
      double pk = h->GetBinCenter(h->GetMaximumBin());
      h->GetXaxis()->SetRange(0, 0);
      double mx = h->GetMaximum();
      TF1 dg("dg", "[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x", pk - 1.6, pk + 1.4);
      dg.SetParameters(pk - 0.3, 0.30, 0.6 * mx, 0.6 * mx, 0.1 * mx, 0.0);
      dg.SetParLimits(0, pk - 1.0, pk + 0.5);
      dg.SetParLimits(1, 0.12, 0.55);
      dg.SetParLimits(2, 0, 2 * mx);
      dg.SetParLimits(3, 0, 2 * mx);
      h->Fit(&dg, "QRN");
      double gs = dg.GetParameter(0), fwhm = 2.3548 * dg.GetParameter(1);
      printf("  %6.1f  %+6.3f   %6.3f\n", E, gs, fwhm);
      gF->SetPoint(ie, E, fwhm);
      if (fwhm < bestF) {
         bestF = fwhm;
         bestE = E;
      }
      delete h;
   }
   printf("\n=== BEST doublet resolution: Ebeam = %.1f MeV, FWHM = %.3f MeV ===\n", bestE, bestF);
   TCanvas *c = new TCanvas("c", "escan2", 800, 600);
   gF->SetTitle("g.s./0.74 doublet FWHM vs beam energy;E_{beam} [MeV];doublet FWHM [MeV]");
   gF->SetMarkerStyle(20);
   gF->Draw("ALP");
   c->SaveAs("pd/plots/escan2.png");
   printf("saved pd/plots/escan2.png\n");
}
