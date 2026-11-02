/// @file plot_pid_range_sim_dp.C
/// @brief RANGE-vs-MOMENTUM PID test on the d + t sims (2026-11-01). For tracks that stop in the gas (range_dp.C cont),
///        the range at fixed momentum scales ~ p^4/m^3, so a triton's is ~0.30 of a deuteron's (ln gap ~1.2), and it does
///        not depend on gain or saturation. Per circle-Brho (bA) slice, d sim vs t sim:
///          dE/dx  : x1 = ln(hs * bA)                (the adopted variable)
///          range  : x2 = ln(range)                  (range = len, AtHelix3D backbone length)
///          both   : Fisher discriminant of (x1, x2), weights from the slice's d and t covariances
///        S = (median_t - median_d) / (hw_t + hw_d) and the misID fractions at the width-weighted midpoint.
///        Containment also requires zfar < zFarMax (955 mm): tracks leaving the drift window end near z 964.
///        Panels: range vs Brho plane (d blue, t red), fraction contained, separation vs Brho, misID vs Brho.
///   root -l -b -q 'plot_pid_range_sim_dp.C("plots/pid_range_sim.png")'
void plot_pid_range_sim_dp(TString png, TString rv = "len", TString R = "/mnt/f/a1975/pid_range/", double zFarMax = 955)
{
   gStyle->SetOptStat(0);
   TChain d("w"), t("w"), dr("r"), trr("r");
   for (int s = 5001; s <= 5004; ++s) { d.Add(Form("/mnt/f/a1975_C16_dd_sim/pv_el_s%d.root", s)); dr.Add(Form("%sr_el_s%d.root", R.Data(), s)); }
   for (const char *s : {"gs", "ex2"}) { t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); trr.Add(Form("%sr_%s.root", R.Data(), s)); }
   d.AddFriend(&dr, "r"); t.AddFriend(&trr, "r");
   TTree *T[2] = {&d, &t};
   TString base = "va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&bA>0&&r.ok";
   // contained = far end inside the gas; range_dp.C's z < 995 lets through tracks that LEAVE the drift window at
   // z ~ 964 mm (sim: 139 of 875 "contained" tritons at bA 0.3-0.6 had len < 300 mm, all with zfar 963-965)
   TString sel = base + Form("&&r.cont&&r.zfar<%g&&r.zfar>10&&r.%s>20", zFarMax, rv.Data());
   auto *c = new TCanvas("c", "c", 1900, 1250);
   c->Divide(2, 2);
   // (1) plane
   c->cd(1); gPad->SetGrid(); gPad->SetLogy();
   auto *fr = gPad->DrawFrame(0.1, 20, 0.9, 3000);
   fr->SetTitle(Form("stopping tracks: range vs circle B#rho;B#rho_{circle} [T m];range (%s) [mm]", rv.Data()));
   const int col[2] = {kBlue + 1, kRed + 1};
   for (int s = 0; s < 2; ++s) {
      long n = T[s]->Draw(Form("r.%s:bA", rv.Data()), sel, "goff");
      auto *g = new TGraph(n, T[s]->GetV2(), T[s]->GetV1());
      g->SetMarkerStyle(1); g->SetMarkerColorAlpha(col[s], 0.6); g->Draw("P");
   }
   const double bl[] = {0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.75};
   const int nb = 7;
   TGraph *gC[2], *gS[3], *gM[3][2];
   for (int s = 0; s < 2; ++s) gC[s] = new TGraph();
   for (int v = 0; v < 3; ++v) { gS[v] = new TGraph(); gM[v][0] = new TGraph(); gM[v][1] = new TGraph(); }
   printf("%-10s %-14s %-24s %-24s %-24s\n", "bA slice", "contained d/t", "dE/dx S t>d d>t", "range S t>d d>t", "both S t>d d>t");
   for (int k = 0; k < nb; ++k) {
      const double x = 0.5 * (bl[k] + bl[k + 1]);
      TString sl = Form("&&bA>%g&&bA<%g", bl[k], bl[k + 1]);
      double fc[2];
      std::vector<std::array<double, 2>> v[2];
      for (int s = 0; s < 2; ++s) {
         double nall = T[s]->GetEntries(base + sl), ncon = T[s]->GetEntries(sel + sl);
         fc[s] = nall > 0 ? ncon / nall : 0;
         gC[s]->SetPoint(gC[s]->GetN(), x, 100 * fc[s]);
         long n = T[s]->Draw(Form("log(hs*bA):log(r.%s)", rv.Data()), sel + sl, "goff");
         for (long i = 0; i < n; ++i) v[s].push_back({T[s]->GetV1()[i], T[s]->GetV2()[i]});
      }
      if (v[0].size() < 15 || v[1].size() < 15) { printf("%.2f-%.2f  N d %zu t %zu\n", bl[k], bl[k + 1], v[0].size(), v[1].size()); continue; }
      // Fisher weights
      double mu[2][2] = {{0, 0}, {0, 0}}, S[2][2] = {{0, 0}, {0, 0}};
      for (int s = 0; s < 2; ++s) {
         for (auto &p : v[s]) { mu[s][0] += p[0]; mu[s][1] += p[1]; }
         mu[s][0] /= v[s].size(); mu[s][1] /= v[s].size();
         for (auto &p : v[s]) {
            double a = p[0] - mu[s][0], b = p[1] - mu[s][1];
            S[0][0] += a * a / v[s].size(); S[0][1] += a * b / v[s].size(); S[1][1] += b * b / v[s].size();
         }
      }
      S[1][0] = S[0][1];
      double det = S[0][0] * S[1][1] - S[0][1] * S[0][1], dm0 = mu[1][0] - mu[0][0], dm1 = mu[1][1] - mu[0][1];
      double w0 = (S[1][1] * dm0 - S[0][1] * dm1) / det, w1 = (-S[1][0] * dm0 + S[0][0] * dm1) / det;
      printf("%.2f-%.2f  %4.0f%% / %4.0f%%  ", bl[k], bl[k + 1], 100 * fc[0], 100 * fc[1]);
      for (int vv = 0; vv < 3; ++vv) {
         std::vector<double> z[2];
         for (int s = 0; s < 2; ++s)
            for (auto &p : v[s]) z[s].push_back(vv == 0 ? p[0] : vv == 1 ? p[1] : w0 * p[0] + w1 * p[1]);
         // range DEcreases with mass: flip so the triton side is always the high side
         if (vv == 1) for (int s = 0; s < 2; ++s) for (auto &q : z[s]) q = -q;
         for (int s = 0; s < 2; ++s) std::sort(z[s].begin(), z[s].end());
         auto q = [](std::vector<double> &a, double f) { return a[std::min(a.size() - 1, (size_t)(f * a.size()))]; };
         double md = q(z[0], .5), mt = q(z[1], .5), hd = .5 * (q(z[0], .84) - q(z[0], .16)), ht = .5 * (q(z[1], .84) - q(z[1], .16));
         double cut = md + (mt - md) * hd / (hd + ht);
         double fd = double(z[0].end() - std::upper_bound(z[0].begin(), z[0].end(), cut)) / z[0].size();
         double ft = double(std::upper_bound(z[1].begin(), z[1].end(), cut) - z[1].begin()) / z[1].size();
         double Sv = (mt - md) / (hd + ht);
         gS[vv]->SetPoint(gS[vv]->GetN(), x, Sv);
         gM[vv][0]->SetPoint(gM[vv][0]->GetN(), x, 100 * ft); gM[vv][1]->SetPoint(gM[vv][1]->GetN(), x, 100 * fd);
         printf("S %4.2f %5.1f%% %5.2f%%     ", Sv, 100 * ft, 100 * fd);
      }
      printf("\n");
   }
   c->cd(2); gPad->SetGrid();
   auto *f2 = gPad->DrawFrame(0.1, 0, 0.8, 105);
   f2->SetTitle("fraction of tracks that stop in the gas;B#rho_{circle} [T m];contained [%]");
   for (int s = 0; s < 2; ++s) { gC[s]->SetMarkerStyle(20); gC[s]->SetMarkerColor(col[s]); gC[s]->SetLineColor(col[s]); gC[s]->SetLineWidth(2); gC[s]->Draw("LP"); }
   auto *l2 = new TLegend(0.6, 0.75, 0.88, 0.88); l2->AddEntry(gC[0], "sim deuterons", "lp"); l2->AddEntry(gC[1], "sim tritons", "lp"); l2->Draw();
   const char *nm[3] = {"dE/dx (adopted)", "range", "dE/dx + range"};
   const int vc[3] = {kBlue + 2, kOrange + 7, kGreen + 2};
   c->cd(3); gPad->SetGrid();
   auto *f3 = gPad->DrawFrame(0.1, 0, 0.8, 8);
   f3->SetTitle("d/t separation, stopping tracks (sim);B#rho_{circle} [T m];#Delta median / sum of half-widths");
   auto *l3 = new TLegend(0.55, 0.72, 0.88, 0.88);
   for (int v = 0; v < 3; ++v) { gS[v]->SetMarkerStyle(20 + v); gS[v]->SetMarkerSize(1.4); gS[v]->SetMarkerColor(vc[v]); gS[v]->SetLineColor(vc[v]); gS[v]->SetLineWidth(2); gS[v]->Draw("LP"); l3->AddEntry(gS[v], nm[v], "lp"); }
   l3->Draw();
   c->cd(4); gPad->SetGrid(); gPad->SetLogy();
   auto *f4 = gPad->DrawFrame(0.1, 0.01, 0.8, 50);
   f4->SetTitle("deuterons in the triton gate (solid), tritons lost (dashed);B#rho_{circle} [T m];misID [%]");
   for (int v = 0; v < 3; ++v)
      for (int s = 0; s < 2; ++s) {
         auto *g = gM[v][s];
         // log scale: floor zero at 0.02 % for display
         for (int i = 0; i < g->GetN(); ++i) if (g->GetY()[i] < 0.02) g->SetPoint(i, g->GetX()[i], 0.02);
         g->SetMarkerStyle(20 + v); g->SetMarkerColor(vc[v]); g->SetLineColor(vc[v]); g->SetLineWidth(2); g->SetLineStyle(s ? 1 : 2);
         g->Draw("LP");
      }
   c->SaveAs(png);
}
