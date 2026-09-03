/// @file gain_report_C15d.C
/// @brief Does the per-run gain match actually work? Four tests, on the data, before any gate.
///
///   root -b -q 'gain_report_C15d.C()'
///
/// The user's requirement is explicit: no gate until the gain matching is shown to work. A gate is
/// a polygon in (sqrt(dE/dx), Brho) drawn ONCE and applied to every run, so if the gain drifts and
/// is not corrected, the same polygon selects a different part of the band in run 13 than in run
/// 133 -- silently, and with a perfectly plausible yield in both.
///
/// TEST 1  SPREAD.      How far does the gain actually drift over the run set? If it is a few
///                      percent the correction hardly matters; if it is tens of percent it does.
/// TEST 2  SMOOTHNESS.  Real gain drift is smooth in time; estimator noise is not. Compare the
///                      mean |change| between CONSECUTIVE runs against the total drift. Noise
///                      comparable to the drift means the factors are measuring nothing.
/// TEST 3  ALIGNMENT.   The point of the correction: does the per-run band position line up after
///                      matching? Scatter of the per-run median across runs, before vs after.
/// TEST 4  RESOLUTION.  Does the pooled band get NARROWER? Alignment can be faked by a factor that
///                      merely rescales; a genuine match also sharpens the summed band. Measured
///                      in single-band Brho slices, because a window straddling two bands measures
///                      their ratio and looks smooth while being wrong.

#include "gain_C15d.h"

namespace {
/// median of a vector, sorted in place
double med(std::vector<double> &v)
{
   if (v.empty()) return 0;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}
/// robust width: half the 16-84 percentile interval, i.e. a Gaussian sigma without the tails
double sig(std::vector<double> &v)
{
   if (v.size() < 20) return 0;
   std::sort(v.begin(), v.end());
   return 0.5 * (v[(size_t)(0.84 * v.size())] - v[(size_t)(0.16 * v.size())]);
}
} // namespace

void gain_report_C15d(TString inDir = "/home/yassid/C15d_reco/", TString gainCsv = "gainmatch_C15d.csv",
                      TString outDir = "plots/", Double_t brLo = 0.30, Double_t brHi = 0.40)
{
   gSystem->Load("libAtTools.so");
   gSystem->mkdir(outDir, kTRUE);
   gStyle->SetOptStat(0);

   auto gain = LoadGainTable_C15d(gainCsv.Data());
   if (gain.empty()) {
      std::cout << "\033[1;31mERROR: no gain table at " << gainCsv << "\033[0m\n";
      return;
   }

   // ---- read every per-run PID cache once -------------------------------------------------
   TString ls = gSystem->GetFromPipe("ls -1 " + inDir + "*_pid.root 2>/dev/null");
   TObjArray *files = ls.Tokenize("\n");
   std::map<int, std::vector<double>> rawAnchor;              // run -> sqrt(dEdx) in the anchor window
   std::map<int, std::vector<std::pair<double, double>>> all; // run -> (brho, sqrt(dEdx)) everywhere
   for (int i = 0; i < files->GetEntries(); ++i) {
      TString fn = ((TObjString *)files->At(i))->GetString();
      TString base = gSystem->BaseName(fn.Data());
      // "run_0013_pid.root": the digits start at index 4, not 3. base(3,4) is "_001", whose
      // Atoi is 0, so EVERY run mapped to run 0 -- the report then said "runs with data : 1"
      // and every test degenerated to a single point with factor 1.0.
      int run = TString(base(4, 4)).Atoi();
      TFile f(fn);
      auto *t = (TTree *)f.Get("pid");
      if (!t) continue;
      Double_t dedx = 0, brho = 0;
      Int_t valid = 1;
      t->SetBranchAddress("dEdx", &dedx);
      t->SetBranchAddress("brho", &brho);
      if (t->GetBranch("valid")) t->SetBranchAddress("valid", &valid);
      for (Long64_t k = 0; k < t->GetEntries(); ++k) {
         t->GetEntry(k);
         if (!valid || dedx <= 0) continue;
         const double s = std::sqrt(dedx);
         all[run].emplace_back(brho, s);
         if (brho >= brLo && brho <= brHi) rawAnchor[run].push_back(s);
      }
   }
   if (rawAnchor.empty()) {
      std::cout << "\033[1;31mERROR: no _pid.root caches in " << inDir << "\033[0m\n";
      return;
   }

   printf("\n\033[1;33m======================  GAIN QUALITY REPORT  ======================\033[0m\n");
   printf("  runs with data : %zu      anchor window Brho [%.2f, %.2f]\n", rawAnchor.size(), brLo, brHi);

   // ---- TEST 1: spread ---------------------------------------------------------------------
   std::vector<int> runs;
   for (auto &kv : rawAnchor) runs.push_back(kv.first);
   std::sort(runs.begin(), runs.end());
   std::vector<double> facs;
   int nMissing = 0;
   for (int r : runs) {
      bool miss = false;
      facs.push_back(GainFactor_C15d(gain, r, miss));
      if (miss) ++nMissing;
   }
   if (nMissing)
      printf("  \033[1;33m  %d of %zu runs have NO gain entry (interpolated or held)\033[0m\n", nMissing,
             runs.size());
   double fmin = *std::min_element(facs.begin(), facs.end());
   double fmax = *std::max_element(facs.begin(), facs.end());
   printf("\n\033[1;36m  TEST 1  SPREAD\033[0m\n");
   printf("    factor %.4f .. %.4f  -> %.1f %% drift across the run set\n", fmin, fmax, 100 * (fmax / fmin - 1));

   // ---- TEST 2: smoothness -----------------------------------------------------------------
   double sumStep = 0, maxStep = 0;
   int nStep = 0;
   for (size_t i = 1; i < runs.size(); ++i) {
      double d = std::fabs(facs[i] - facs[i - 1]);
      sumStep += d;
      maxStep = std::max(maxStep, d);
      ++nStep;
   }
   const double meanStep = nStep ? sumStep / nStep : 0;
   const double total = fmax - fmin;
   printf("\n\033[1;36m  TEST 2  SMOOTHNESS\033[0m\n");
   printf("    mean |change| between consecutive runs : %.4f\n", meanStep);
   printf("    max  |change|                          : %.4f\n", maxStep);
   printf("    total drift                            : %.4f\n", total);
   printf("    step / drift = %.3f  -> %s\n", total > 0 ? meanStep / total : 0,
          (total > 0 && meanStep / total < 0.15)
             ? "\033[1;32mSMOOTH: real gain drift, not estimator noise\033[0m"
             : "\033[1;31mNOISY: the step-to-step scatter rivals the drift; the factors may be noise\033[0m");

   // ---- TEST 3: alignment ------------------------------------------------------------------
   std::vector<double> medRaw, medCor;
   for (int r : runs) {
      auto v = rawAnchor[r];
      if (v.size() < 100) continue;
      bool miss = false;
      const double f = GainFactor_C15d(gain, r, miss);
      std::vector<double> w;
      w.reserve(v.size());
      for (double s : v) w.push_back(s * std::sqrt(f));
      medRaw.push_back(med(v));
      medCor.push_back(med(w));
   }
   auto rms = [](std::vector<double> v) {
      if (v.size() < 2) return 0.0;
      double m = 0;
      for (double x : v) m += x;
      m /= v.size();
      double s = 0;
      for (double x : v) s += (x - m) * (x - m);
      return std::sqrt(s / v.size()) / m * 100.0; // % of the mean
   };
   const double aRaw = rms(medRaw), aCor = rms(medCor);
   printf("\n\033[1;36m  TEST 3  ALIGNMENT  (scatter of the per-run band position, %zu runs)\033[0m\n", medRaw.size());
   printf("    raw       : %.2f %%\n", aRaw);
   printf("    matched   : %.2f %%   -> %s\n", aCor,
          aCor < aRaw ? Form("\033[1;32m%.1fx tighter\033[0m", aRaw / std::max(1e-9, aCor))
                      : "\033[1;31mNO IMPROVEMENT\033[0m");

   // ---- TEST 4: resolution, in several single-band slices ----------------------------------
   printf("\n\033[1;36m  TEST 4  RESOLUTION  (pooled band width, sigma/median)\033[0m\n");
   printf("    %-16s %10s %10s %10s\n", "Brho slice", "raw", "matched", "change");
   const double slices[][2] = {{0.25, 0.30}, {0.30, 0.40}, {0.40, 0.50}, {0.50, 0.60}};
   for (auto &s : slices) {
      std::vector<double> R, C;
      for (int r : runs) {
         bool miss = false;
         const double f = GainFactor_C15d(gain, r, miss);
         for (auto &pr : all[r])
            if (pr.first >= s[0] && pr.first < s[1]) {
               R.push_back(pr.second);
               C.push_back(pr.second * std::sqrt(f));
            }
      }
      if (R.size() < 500) continue;
      const double wR = sig(R) / std::max(1e-9, med(R)), wC = sig(C) / std::max(1e-9, med(C));
      printf("    %-16s %9.2f%% %9.2f%% %9.1f%%%s\n", Form("%.2f-%.2f", s[0], s[1]), 100 * wR, 100 * wC,
             100 * (wC / wR - 1), wC < wR ? "  narrower" : "  WIDER");
   }

   // ---- plots -------------------------------------------------------------------------------
   auto *gF = new TGraph(), *gR = new TGraph(), *gC = new TGraph();
   for (size_t i = 0, k = 0; i < runs.size(); ++i) {
      gF->SetPoint(gF->GetN(), runs[i], facs[i]);
      if (rawAnchor[runs[i]].size() >= 100 && k < medRaw.size()) {
         gR->SetPoint(gR->GetN(), runs[i], medRaw[k]);
         gC->SetPoint(gC->GetN(), runs[i], medCor[k]);
         ++k;
      }
   }
   auto *c = new TCanvas("cg", "gain", 1500, 560);
   c->Divide(2, 1);
   c->cd(1);
   gF->SetTitle("per-run gain factor;run;factor");
   gF->SetMarkerStyle(20);
   gF->SetLineWidth(2);
   gF->Draw("ALP");
   c->cd(2);
   gR->SetTitle(Form("band position in Brho %.2f-%.2f;run;median #sqrt{dE/dx}", brLo, brHi));
   gR->SetMarkerStyle(20);
   gR->SetMarkerColor(kGray + 2);
   gR->SetLineColor(kGray + 2);
   gR->Draw("ALP");
   gC->SetMarkerStyle(20);
   gC->SetMarkerColor(kAzure + 2);
   gC->SetLineColor(kAzure + 2);
   gC->SetLineWidth(2);
   gC->Draw("LP same");
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.04);
   tx->SetTextColor(kGray + 2);
   tx->DrawLatex(0.45, 0.86, Form("raw  (scatter %.2f%%)", aRaw));
   tx->SetTextColor(kAzure + 2);
   tx->DrawLatex(0.45, 0.81, Form("matched  (scatter %.2f%%)", aCor));
   c->SaveAs(outDir + "gain_report_C15d.png");

   printf("\n  wrote %sgain_report_C15d.png\n", outDir.Data());
   printf("\n\033[1;33m  VERDICT: %s\033[0m\n",
          (aCor < aRaw && total > 0 && meanStep / total < 0.15)
             ? "the gain match works -- the plane is safe to draw a gate on"
             : "DO NOT GATE YET -- see the failing test above");
   printf("\033[1;33m==================================================================\033[0m\n\n");
}
