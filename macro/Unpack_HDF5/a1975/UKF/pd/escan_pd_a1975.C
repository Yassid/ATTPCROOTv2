/// @file escan_pd_a1975.C
/// @brief Beam-energy scan for 16C(p,d)15C: find the Ebeam that best resolves the
///        g.s./0.74 doublet and flattens the Ex-vs-theta_cm dependence.
///
/// A wrong effective beam energy introduces a theta_cm-dependent Ex error that
/// smears the peaks. For each Ebeam (default 189-195 MeV, 100 keV steps) this
/// recomputes Ex (and theta_cm) per track from the cached (ke, theta), fits the
/// g.s./0.74 doublet with the separation FIXED at 0.740 MeV (stable), and reports:
///   - g.s. peak FWHM  (resolution; smaller = better)
///   - the band spread across theta_cm (flatness; smaller = flatter)
///   - a free-separation fit (to check the doublet gap ~ 0.740)
/// Picks the Ebeam minimizing the g.s. FWHM.
///
///   root -b -q 'pd/escan_pd_a1975.C'
///   root -b -q 'pd/escan_pd_a1975.C(190.0,194.0,0.1)'

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static const double U = 931.49401, M1 = 16.0147 * U, M2 = 1.007825 * U, M3 = 2.01410178 * U, M4 = 15.0105993 * U;
static void kine(double Kp, double thR, double Ke, double &Ex, double &thcm)
{
   double Et1 = Kp + M1, Et3 = Ke + M3, Et4 = Et1 + M2 - Et3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double m4e = std::sqrt((std::cos(thR) * omega2(s, M1 * M1, M2 * M2) * omega2(uu, M2 * M2, M3 * M3) -
                           (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                            (2 * M2 * M2) +
                          s + uu - M2 * M2);
   Ex = m4e - M4;
   double tt = M2 * M2 + m4e * m4e - 2 * M2 * Et4;
   thcm = (TMath::Pi() -
           std::acos((s * s + s * (2 * tt - M1 * M1 - M2 * M2 - M3 * M3 - m4e * m4e) + (M1 * M1 - M2 * M2) *
                                                                                          (M3 * M3 - m4e * m4e)) /
                     (omega2(s, M1 * M1, M2 * M2) * omega2(s, M3 * M3, m4e)))) *
          TMath::RadToDeg();
}

void escan_pd_a1975(double Elo = 189.0, double Ehi = 195.0, double dstep = 0.1, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ke, th, ex, thcm;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("thcm", &thcm);
   std::vector<float> vke, vth, vtc;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18 && ke > 0 && ke < 200 && th > 0 && th < 90) {
         vke.push_back(ke);
         vth.push_back(th);
         vtc.push_back(thcm); // cached theta_cm (binning axis for flatness)
      }
   }
   f->Close();
   printf("scan tracks: %zu\n\n", vke.size());

   int nE = (int)std::round((Ehi - Elo) / dstep) + 1;
   TGraph *gF = new TGraph(), *gS = new TGraph(), *gSep = new TGraph();
   printf("  Ebeam   gs_pos   gs_FWHM   sep(free)   band_spread(thcm)\n");
   double bestF = 1e9, bestE = Elo;
   for (int ie = 0; ie < nE; ++ie) {
      double E = Elo + ie * dstep;
      TH1F *h = new TH1F("h", "", 200, -6, 14);
      h->SetDirectory(nullptr);
      TProfile *pr = new TProfile("pr", "", 16, 20, 140, -3, 3);
      pr->SetDirectory(nullptr);
      for (size_t i = 0; i < vke.size(); ++i) {
         double Ex, tc;
         kine(E, vth[i] * TMath::DegToRad(), vke[i], Ex, tc);
         if (std::isnan(Ex))
            continue;
         h->Fill(Ex);
         if (Ex > -1.0 && Ex < 1.6)
            pr->Fill(vtc[i], Ex); // g.s.-region band vs CACHED theta_cm (flatness)
      }
      // locate doublet: tallest peak in [-2,2.5] is the 0.74 peak; g.s. ~ that-0.74
      h->GetXaxis()->SetRangeUser(-2, 2.5);
      double pk = h->GetBinCenter(h->GetMaximumBin());
      h->GetXaxis()->SetRange(0, 0);
      // fixed-0.74-separation doublet fit
      TF1 dg("dg", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1]-0.740)/[4])^2)", pk - 1.1, pk + 1.1);
      double mx = h->GetMaximum();
      dg.SetParameters(0.6 * mx, pk - 0.4, 0.30, mx, 0.30);
      dg.SetParLimits(1, pk - 1.1, pk + 0.3);
      dg.SetParLimits(2, 0.10, 1.0);
      dg.SetParLimits(4, 0.10, 1.0);
      h->Fit(&dg, "QRN");
      double gs = dg.GetParameter(1), fwhm = 2.3548 * dg.GetParameter(2);
      // free-separation fit for the gap
      TF1 fg("fg", "gaus(0)+gaus(3)", gs - 0.9, gs + 1.7);
      fg.SetParameters(dg.GetParameter(0), gs, 0.3, dg.GetParameter(3), gs + 0.74, 0.3);
      fg.SetParLimits(1, gs - 0.6, gs + 0.4);
      fg.SetParLimits(4, gs + 0.3, gs + 1.3);
      h->Fit(&fg, "QRN");
      double sep = fabs(fg.GetParameter(4) - fg.GetParameter(1));
      // band spread across theta_cm (stddev of profile bin means, populated bins)
      std::vector<double> mns;
      for (int b = 1; b <= pr->GetNbinsX(); ++b)
         if (pr->GetBinEntries(b) > 20)
            mns.push_back(pr->GetBinContent(b));
      double mean = 0;
      for (double m : mns)
         mean += m;
      mean /= std::max((size_t)1, mns.size());
      double sd = 0;
      for (double m : mns)
         sd += (m - mean) * (m - mean);
      sd = std::sqrt(sd / std::max((size_t)1, mns.size()));

      printf("  %6.1f  %+6.3f   %6.3f    %6.3f      %6.3f\n", E, gs, fwhm, sep, sd);
      gF->SetPoint(ie, E, fwhm);
      gS->SetPoint(ie, E, sd);
      gSep->SetPoint(ie, E, sep);
      if (fwhm < bestF) {
         bestF = fwhm;
         bestE = E;
      }
      delete h;
      delete pr;
   }
   printf("\n=== BEST (min g.s. FWHM): Ebeam = %.1f MeV, FWHM = %.3f MeV ===\n", bestE, bestF);

   TCanvas *c = new TCanvas("c", "escan", 1200, 500);
   c->Divide(3, 1);
   c->cd(1);
   gF->SetTitle("g.s. FWHM vs Ebeam;Ebeam [MeV];g.s. FWHM [MeV]");
   gF->SetMarkerStyle(20);
   gF->Draw("ALP");
   c->cd(2);
   gS->SetTitle("Ex band spread vs theta_{cm};Ebeam [MeV];band spread [MeV]");
   gS->SetMarkerStyle(20);
   gS->Draw("ALP");
   c->cd(3);
   gSep->SetTitle("doublet separation;Ebeam [MeV];sep [MeV]");
   gSep->SetMarkerStyle(20);
   gSep->Draw("ALP");
   TLine *l074 = new TLine(Elo, 0.740, Ehi, 0.740);
   l074->SetLineColor(kRed);
   l074->SetLineStyle(2);
   l074->Draw();
   c->SaveAs("pd/plots/escan.png");
   printf("saved pd/plots/escan.png\n");
}
