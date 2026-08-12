/// @file xsec_pd.C
/// @brief Absolute 16C(p,d)15C ground-state cross-section: yield / acceptance / luminosity.
///
///     dsigma/dOmega [mb/sr] = N(theta_cm) / dOmega / acceptance(theta_cm) / L
///
/// L = 316.4 mb^-1, measured from the 16C(p,p) elastic against DWBA (pp/dwba_compare_pp.C). It
/// transfers here unchanged: same runs, same IC gate, same protons in the same gas over the same
/// full ~940 mm of target, so both the beam count and the areal density are identical. It is a
/// LUMINOSITY, not a fudge factor -- 316.4 mb^-1 = 3.16e29 cm^-2 = 1.70e8 gated 16C ions on
/// 1.86e21 protons/cm2.
///
/// ACCEPTANCE is per state, from macro/Simulation/ATTPC/16C_pd/diagnostics/acceptance_<tag>.root.
/// It has to be per state because each Q-value puts the deuteron at a different energy for the
/// same theta_cm. It is a reco+fit acceptance with NO PID efficiency folded in, which is correct
/// to about 1 %: a gate drawn on the simulated plane keeps 99.3 % of truth-matched deuterons, so
/// the deuteron PID selection is essentially lossless. (Applying the DATA gate to the simulation
/// gives 79.8 %, but that measures the offset between the two planes, not an efficiency.)
///
/// BINS BELOW accMin ARE DROPPED, NOT CORRECTED, and the acceptance is required to be at least as
/// fine as the yield -- dividing a yield by a coarser bin's average once corrected a bin whose
/// true acceptance was zero.
///
/// The g.s. is a PEAK here, not the whole spectrum, so the yield is the count inside an Ex window.
/// Where neighbouring states leak into that window the point is an upper limit; the printed table
/// gives the raw count per bin so thin bins are visible rather than implied.
///
///   root -b -q 'pp/xsec_pd.C("/path/pd_kin.root")'

static double xd_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
/// identical chain and masses to acceptance_C16pd.C, so yield and acceptance share an axis
static std::pair<double, double> xd_kine(double Kp, double thl, double Ke)
{
   const double u = 931.49401;
   double m1 = 16.0147 * u, m2 = 1.007825 * u, m3 = 2.014102 * u, m4 = 15.010599 * u;
   double Et1 = Kp + m1, Et3 = Ke + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * xd_om2(s, m1 * m1, m2 * m2) * xd_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   if (a < 0) return {std::nan(""), std::nan("")};
   double m4x = std::sqrt(a), ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * Et4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
               (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (xd_om2(s, m1 * m1, m2 * m2) * xd_om2(s, m3 * m3, m4x * m4x));
   if (c < -1 || c > 1) return {std::nan(""), std::nan("")};
   return {ex, (TMath::Pi() - std::acos(c)) * TMath::RadToDeg()};
}

void xsec_pd(TString cache, TString accFile = "../../../../Simulation/ATTPC/16C_pd/diagnostics/acceptance_gs.root",
             TString accHist = "hAcc_gs", Double_t lumi = 316.4, Double_t exLo = -0.6, Double_t exHi = 0.6,
             Double_t exShift = -0.38, Double_t Ebeam = 192.0, Double_t chi2max = 5.0, Double_t icMin = 1000,
             Double_t icMax = 1350, Double_t accMin = 0.15, TString tag = "gs")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fa = TFile::Open(here + "/" + accFile);
   TH1D *A = fa && !fa->IsZombie() ? (TH1D *)fa->Get(accHist) : nullptr;
   if (!A) { printf("\033[1;31mno '%s' in %s\033[0m\n", accHist.Data(), accFile.Data()); return; }
   A = (TH1D *)A->Clone("Acl"); A->SetDirectory(nullptr); fa->Close();
   const double accW = A->GetBinWidth(1);
   printf("\n  acceptance %s: %d bins of %.1f deg\n", accHist.Data(), A->GetNbinsX(), accW);

   TFile *f = TFile::Open(cache);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("pk") : nullptr;
   if (!t) { printf("\033[1;31mcannot open %s\033[0m\n", cache.Data()); return; }
   float ke, th, vz, c2, ic;
   t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vz", &vz);
   t->SetBranchAddress("chi2ndf", &c2); t->SetBranchAddress("ic", &ic);

   // yield binned exactly like the acceptance -- never coarser
   const int NB = A->GetNbinsX();
   auto *hY = new TH1D("hY", "", NB, A->GetXaxis()->GetXmin(), A->GetXaxis()->GetXmax());
   hY->Sumw2();
   // Detach from the cache file: histograms created while a TFile is open belong to that
   // directory, and f->Close() below would delete them out from under the code that follows.
   hY->SetDirectory(nullptr);
   auto *hEx = new TH1D("hEx", "16C(p,d)^{15}C;E_{x} [MeV];counts", 200, -3, 11);
   hEx->SetDirectory(nullptr);
   auto *h2 = new TH2D("h2", "selected;#theta_{cm} [deg];E_{x} [MeV]", NB, A->GetXaxis()->GetXmin(),
                       A->GetXaxis()->GetXmax(), 160, -3, 11);
   h2->SetDirectory(nullptr);
   long nAll = 0, nSel = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2max && ic > icMin && ic < icMax)) continue;
      auto [ex0, tcm] = xd_kine(Ebeam, th * TMath::DegToRad(), ke);
      if (std::isnan(ex0)) continue;
      double ex = ex0 + exShift;
      ++nAll; hEx->Fill(ex); h2->Fill(tcm, ex);
      if (ex < exLo || ex > exHi) continue;
      ++nSel; hY->Fill(tcm);
   }
   f->Close();
   printf("  %ld deuterons pass chi2/IC, %ld inside Ex [%.2f, %.2f] (%.1f %%)\n\n", nAll, nSel, exLo, exHi,
          100.0 * nSel / std::max(1L, nAll));

   auto dOm = [](double lo, double hi) {
      return 2 * TMath::Pi() * (std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad()));
   };
   auto *hX = (TH1D *)hY->Clone("hX"); hX->Reset(); hX->SetDirectory(nullptr);

   printf("  theta_cm |  counts |   acc  +- err |  dsigma/dOmega [mb/sr]\n");
   int nDrop = 0;
   for (int b = 1; b <= NB; ++b) {
      double lo = hY->GetBinLowEdge(b), hi = lo + hY->GetBinWidth(b), c = hY->GetBinCenter(b);
      double n = hY->GetBinContent(b), en = hY->GetBinError(b);
      if (n <= 0) continue;
      int ab = A->FindBin(c);
      double a = A->GetBinContent(ab), ea = A->GetBinError(ab);
      double dO = dOm(lo, hi);
      if (a < accMin) { printf("  %3.0f-%3.0f  | %7.0f | %.3f       |  DROPPED (acc < %.2f)\n", lo, hi, n, a, accMin);
                        ++nDrop; continue; }
      double v = n / dO / a / lumi;
      double e = v * std::sqrt(std::pow(en / n, 2) + std::pow(ea / std::max(a, 1e-9), 2));
      hX->SetBinContent(b, v); hX->SetBinError(b, e);
      printf("  %3.0f-%3.0f  | %7.0f | %.3f %.3f |  %9.4g +- %.3g\n", lo, hi, n, a, ea, v, e);
   }
   if (nDrop) printf("\n  %d bins dropped for acceptance below %.2f\n", nDrop, accMin);

   TCanvas *c1 = new TCanvas("cx", "pd cross section", 1500, 560);
   c1->Divide(3, 1);
   c1->cd(1); gPad->SetLogy();
   hEx->SetLineWidth(2); hEx->Draw("hist");
   auto *l1 = new TLine(exLo, 0, exLo, hEx->GetMaximum()); l1->SetLineColor(kRed + 1); l1->SetLineStyle(2); l1->Draw();
   auto *l2 = new TLine(exHi, 0, exHi, hEx->GetMaximum()); l2->SetLineColor(kRed + 1); l2->SetLineStyle(2); l2->Draw();
   c1->cd(2); gPad->SetLogz(); h2->Draw("colz");
   c1->cd(3); gPad->SetLogy(); gPad->SetGridy();
   hX->SetTitle(TString::Format("^{16}C(p,d)^{15}C %s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", tag.Data()));
   hX->SetMarkerStyle(20); hX->SetMarkerSize(1.1); hX->SetLineWidth(2);
   hX->GetXaxis()->SetRangeUser(0, 180);
   hX->Draw("E1");

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/xsec_pd_" + tag + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/xsec_pd_" + tag + ".root", "RECREATE");
   hX->Write("dsdo_mb_sr"); hY->Write("yield"); A->Write("acceptance"); hEx->Write("ex"); h2->Write("ex_vs_thcm");
   fo.Close();
   printf("\n  wrote plots/xsec_pd_%s.png/.root\n\n", tag.Data());
}
