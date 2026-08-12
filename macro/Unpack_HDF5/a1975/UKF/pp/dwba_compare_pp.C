/// @file dwba_compare_pp.C
/// @brief Compare the acceptance-corrected 16C(p,p) g.s. angular distribution against DWBA.
///
/// The DWBA file carries five optical-potential variants in columns 2-6 (header labels VV KK MM GG
/// PP) at ELab = 11.50 MeV/u, against a beam of 188 MeV / 16 = 11.75 MeV/u. Column 1 is theta_cm on
/// the SAME convention as the data: dsigma/dOmega diverges as theta_cm -> 0, the Rutherford limit at
/// small momentum transfer, which is the same end as the lowest-energy recoil protons in the data.
///
/// UNITS. The DWBA columns are mb/sr, verified against Rutherford rather than assumed: at 1 deg the
/// pure-Coulomb value for p+16C at Ecm = 10.8 MeV is 6.867e7 mb/sr against the file's 6.83e7, and by
/// 5 deg the file falls below Rutherford as absorption sets in. So the fitted factor CALIBRATES the
/// data: dividing the measured counts/sr by it puts the measurement on an absolute mb/sr axis. The
/// data is plotted divided, the theory in its own units -- scaling the theory UP into arbitrary data
/// units instead would make the axis meaningless.
///
/// THE DATA IS A SHAPE, so one multiplicative factor per curve is fitted rather than assumed. The
/// factor is the error-weighted least-squares scale over fitLo..fitHi only -- the plateau, where the
/// acceptance correction is 1.00-1.03 and trustworthy. Fitting over the full range would let the
/// forward bins, which carry the largest counts and the largest correction, set the normalisation
/// for everything else.
///
/// Agreement is quoted as the rms of ln(data/theory) over the fit window, not chi2/ndf: the data
/// errors here are counting errors on an uncalibrated shape and do not include the acceptance
/// systematic, so a chi2 would look precise for the wrong reason. The rms of the log ratio is
/// scale-free and says how far off the SHAPE is, which is the only thing this comparison can test.
///
///   root -b -q 'pp/dwba_compare_pp.C()'

void dwba_compare_pp(TString dataFile = "plots/angdist_gs_acc_hand5deg.root", TString dataHist = "dsdo_corrected",
                     TString dwbaFile = "/mnt/c/Users/Yassid/Downloads/ATTPCROOTv2-UKF-analysis-12Be_pd/"
                                        "DWBA calculation/16C_pp_0M_DWBA.txt",
                     Double_t fitLo = 20, Double_t fitHi = 120, TString png = "plots/dwba_compare_pp.png")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + dataFile);
   TH1D *D = fd && !fd->IsZombie() ? (TH1D *)fd->Get(dataHist) : nullptr;
   if (!D) { printf("\033[1;31mno '%s' in %s\033[0m\n", dataHist.Data(), dataFile.Data()); return; }
   D = (TH1D *)D->Clone("Dcl"); D->SetDirectory(nullptr); fd->Close();

   // ---- DWBA: angle + five potentials ----
   const int NC = 5;
   const char *lab[NC] = {"VV", "KK", "MM", "GG", "PP"};
   std::vector<double> ang; std::vector<std::vector<double>> sig(NC);
   {
      std::ifstream in(dwbaFile.Data());
      if (!in) { printf("\033[1;31mcannot open %s\033[0m\n", dwbaFile.Data()); return; }
      std::string line;
      while (std::getline(in, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream is(line);
         double a; if (!(is >> a)) continue;              // header lines fail here and are skipped
         std::vector<double> v; double y;
         while (is >> y) v.push_back(y);
         if ((int)v.size() < NC) continue;
         // the 0 deg row is the Coulomb singularity (1e39) and is not a data point
         if (a < 0.5) continue;
         ang.push_back(a);
         for (int c = 0; c < NC; ++c) sig[c].push_back(v[c]);
      }
   }
   if (ang.size() < 10) { printf("\033[1;31monly %zu DWBA rows parsed\033[0m\n", ang.size()); return; }
   printf("\n  DWBA: %zu angles, %.1f-%.1f deg, %d potentials\n", ang.size(), ang.front(), ang.back(), NC);

   auto theoryAt = [&](int c, double a) { // linear interpolation onto the data bin centre
      if (a <= ang.front()) return sig[c].front();
      if (a >= ang.back()) return sig[c].back();
      size_t i = std::lower_bound(ang.begin(), ang.end(), a) - ang.begin();
      double t = (a - ang[i - 1]) / (ang[i] - ang[i - 1]);
      return sig[c][i - 1] * (1 - t) + sig[c][i] * t;
   };

   printf("\n  potential   scale (data/theory)   rms ln(data/theory) over %.0f-%.0f deg   points\n", fitLo, fitHi);
   std::vector<TGraph *> gr(NC);
   double best = 1e18; int bestC = 0;
   std::vector<double> scale(NC, 1.0), rms(NC, 0.0);
   for (int c = 0; c < NC; ++c) {
      // error-weighted least-squares scale: minimise sum w (d - s*t)^2  ->  s = sum(w d t)/sum(w t^2)
      double num = 0, den = 0; int n = 0;
      for (int b = 1; b <= D->GetNbinsX(); ++b) {
         double a = D->GetBinCenter(b), d = D->GetBinContent(b), e = D->GetBinError(b);
         if (d <= 0 || e <= 0 || a < fitLo || a > fitHi) continue;
         double t = theoryAt(c, a); if (t <= 0) continue;
         double w = 1.0 / (e * e);
         num += w * d * t; den += w * t * t; ++n;
      }
      if (den <= 0) continue;
      scale[c] = num / den;
      double s2 = 0; int m = 0;
      for (int b = 1; b <= D->GetNbinsX(); ++b) {
         double a = D->GetBinCenter(b), d = D->GetBinContent(b);
         if (d <= 0 || a < fitLo || a > fitHi) continue;
         double t = scale[c] * theoryAt(c, a); if (t <= 0) continue;
         double l = std::log(d / t); s2 += l * l; ++m;
      }
      rms[c] = m ? std::sqrt(s2 / m) : 0;
      printf("  %-10s  %12.4g        %8.3f                        %3d\n", lab[c], scale[c], rms[c], m);
      if (rms[c] < best) { best = rms[c]; bestC = c; }

      gr[c] = new TGraph();
      for (size_t i = 0; i < ang.size(); ++i)
         if (ang[i] >= 5 && ang[i] <= 165) gr[c]->SetPoint(gr[c]->GetN(), ang[i], sig[c][i]);
   }
   printf("\n  best shape agreement: %s (rms %.3f)\n\n", lab[bestC], rms[bestC]);

   printf("\n  calibration: %.4g counts/sr per mb/sr  (using %s)\n\n", scale[bestC], lab[bestC]);
   printf("  theta_cm   data [mb/sr]   %s [mb/sr]   data/theory\n", lab[bestC]);
   for (int b = 1; b <= D->GetNbinsX(); ++b) {
      double a = D->GetBinCenter(b), d = D->GetBinContent(b);
      if (d <= 0) continue;
      double t = theoryAt(bestC, a), dm = d / scale[bestC];
      printf("  %5.1f   %12.4g   %12.4g   %10.2f%s\n", a, dm, t, t > 0 ? dm / t : 0,
             (a < fitLo || a > fitHi) ? "   (outside fit)" : "");
   }

   // put the measurement on the absolute axis rather than dragging the theory onto an arbitrary one
   for (int b = 1; b <= D->GetNbinsX(); ++b) {
      D->SetBinContent(b, D->GetBinContent(b) / scale[bestC]);
      D->SetBinError(b, D->GetBinError(b) / scale[bestC]);
   }

   TCanvas *c1 = new TCanvas("cdw", "DWBA vs data", 1000, 750);
   gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   D->SetTitle("^{16}C(p,p) g.s. vs DWBA;#theta_{cm} [deg];d#sigma/d#Omega  [mb/sr]");
   D->SetMarkerStyle(20); D->SetMarkerSize(1.2); D->SetLineWidth(2); D->SetLineColor(kBlack);
   D->SetMarkerColor(kBlack);
   D->GetXaxis()->SetRangeUser(0, 160);
   D->Draw("E1");
   int col[NC] = {kAzure + 2, kRed + 1, kGreen + 3, kOrange + 7, kMagenta + 2};
   auto *lg = new TLegend(0.58, 0.62, 0.89, 0.88);
   lg->AddEntry(D, Form("data (corrected, /%.4g)", scale[bestC]), "lp");
   for (int c = 0; c < NC; ++c) {
      if (!gr[c]) continue;
      gr[c]->SetLineColor(col[c]); gr[c]->SetLineWidth(c == bestC ? 4 : 2);
      gr[c]->SetLineStyle(c == bestC ? 1 : 2);
      gr[c]->Draw("L same");
      lg->AddEntry(gr[c], Form("%s  (rms %.2f)", lab[c], rms[c]), "l");
   }
   D->Draw("E1 same");
   lg->Draw();
   gSystem->mkdir(here + "/plots", kTRUE);
   c1->SaveAs(here + "/" + png);
   printf("\n  wrote %s\n\n", png.Data());
}
