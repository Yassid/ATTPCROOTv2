/// @file order_score_find_dp.C
/// @brief Backward (d,p) protons, runs 0016 0026 0058 0088: the adopted find805 chain on production clusters (stored September production, and today's code as a
///        control) and on AtPRA::SetRemoveStrayHits clusters. Adopted selection, Ex as in order_score_p_dp.C. Tracks are read through the PID results, so the
///        repeated entries of the genfit trees are not counted twice. Reports the whole gated sample, the proton-band tracks, and the set selected in both of a pair.
///   root -l -b -q 'order_score_find_dp.C("plots/order_dp_find.png")'
namespace ordf {
struct K { double ke, th, chi, vz, vr; };
using Key = std::tuple<int, Long64_t, int>;
double exDP(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613, mp = 938.272, M4 = mC + md - mp + 1.496;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ep = ke + mp, pp = std::sqrt(Ep * Ep - mp * mp);
   const double E4 = Eb + md - Ep, px = -pp * std::sin(th), pz = pb - pp * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - M4;
}
bool selected(const K &k, double &ex)
{
   const double deg = k.th * 57.29578;
   if (deg < 89 || deg > 161 || k.ke <= 1.2 || k.chi >= 13.5 || k.vz < -100 || k.vz > 940) return false;
   ex = exDP(k.ke + 0.25, k.th, 184.25 - 0.01103 * (k.vz - 500));
   return true;
}
void load(TString dir, const char *run, int irun, std::map<Key, K> &out)
{
   TFile f(dir + "run_" + run + "_multifit_genfitter_p.root");
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te); t->SetBranchAddress("AtPIDEvent", &pe);
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!pe->GetEntries() || !te->GetEntries()) continue;
      auto *pid = (AtPIDEvent *)pe->At(0);
      if (!pid || pid->GetSpyral().empty()) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue;
         auto &kx = ft->GetKinematicsXtr();
         const double ndf = ft->GetTrackMetadata()->GetNdf();
         if (kx.kineticEnergy <= 0 || ndf <= 0) continue;
         auto v = ft->GetVertex();
         out[{irun, i, ft->GetTrackID()}] = {kx.kineticEnergy, kx.theta, ft->GetTrackMetadata()->GetChi2() / ndf, v.Z(), std::hypot(v.X(), v.Y())};
      }
   }
}
void line(const char *nm, const std::map<Key, K> &m, const std::set<Key> *only, TH1D *h)
{
   long fit = 0, sel = 0, bnd = 0, junk = 0, off = 0;
   for (auto &e : m) {
      if (only && !only->count(e.first)) continue;
      ++fit;
      double ex;
      if (!selected(e.second, ex)) continue;
      ++sel; bnd += ex > -1.5 && ex < 1; junk += ex < -2; off += e.second.vr > 20;
      if (h) h->Fill(ex);
   }
   printf("  %-30s fitted %6ld  selected %5ld  bound %4ld  junk %4ld (%.1f %%)  off-axis %4ld (%.1f %%)\n", nm, fit, sel, bnd, junk, 100. * junk / std::max(1L, sel), off, 100. * off / std::max(1L, sel));
}
}
void order_score_find_dp(TString png = "plots/order_dp_find.png", TString Q = "/mnt/f/a1975/pra_test/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   gStyle->SetOptStat(0);
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   const int NP = 3;
   const char *nm[NP] = {"stored find805 (September)", "production+FIND today", "guard+FIND"};
   TString dirs[NP] = {"/mnt/f/a1975/gf_dp_find805/", Q + "gf_ctrl_find/", Q + "gf_guard_find/"};
   std::map<ordf::Key, ordf::K> M[NP];
   std::set<ordf::Key> band;
   for (int r = 0; r < 4; ++r) {
      for (int p = 0; p < NP; ++p) ordf::load(dirs[p], runs[r], r, M[p]);
      std::ifstream il((Q + "pband_" + runs[r] + ".txt").Data());
      Long64_t e; int t;
      while (il >> e >> t) band.insert({r, e, t});
   }
   printf("whole gated sample:\n");
   TH1D *h[NP];
   for (int p = 0; p < NP; ++p) { h[p] = new TH1D(Form("h%d", p), ";E_{x} [MeV];tracks / 0.2 MeV", 70, -5, 9); ordf::line(nm[p], M[p], nullptr, h[p]); }
   printf("proton-band tracks (%zu):\n", band.size());
   for (int p = 0; p < NP; ++p) ordf::line(nm[p], M[p], &band, nullptr);
   // control vs guard, track by track
   long both = 0, onlyC = 0, onlyG = 0, stay = 0, gainB = 0, lossB = 0, junkC = 0, junkG = 0, offC = 0, offG = 0;
   std::vector<double> dEx;
   for (auto &e : M[1]) {
      double xc, xg;
      if (!ordf::selected(e.second, xc)) continue;
      auto it = M[2].find(e.first);
      if (it == M[2].end() || !ordf::selected(it->second, xg)) { ++onlyC; lossB += xc > -1.5 && xc < 1; continue; }
      ++both; junkC += xc < -2; junkG += xg < -2; offC += e.second.vr > 20; offG += it->second.vr > 20;
      const bool bc = xc > -1.5 && xc < 1, bg = xg > -1.5 && xg < 1;
      stay += bc && bg; gainB += !bc && bg; lossB += bc && !bg;
      dEx.push_back(std::abs(xg - xc));
   }
   for (auto &e : M[2]) {
      double xc, xg;
      if (!ordf::selected(e.second, xg)) continue;
      auto it = M[1].find(e.first);
      if (it == M[1].end() || !ordf::selected(it->second, xc)) { ++onlyG; gainB += xg > -1.5 && xg < 1; }
   }
   std::sort(dEx.begin(), dEx.end());
   long moved = std::count_if(dEx.begin(), dEx.end(), [](double d) { return d > 0.5; });
   printf("control vs guard+FIND, track by track: selected in both %ld, control only %ld, guard only %ld\n", both, onlyC, onlyG);
   printf("  bound: kept %ld, gained %ld, lost %ld (entering or leaving the selection included)\n", stay, gainB, lossB);
   printf("  common set: junk %ld -> %ld, off-axis %ld -> %ld, median |dEx| %.3f MeV, |dEx|>0.5 on %ld tracks\n", junkC, junkG, offC, offG, dEx[dEx.size() / 2], moved);
   auto *c = new TCanvas("c", "c", 1100, 700);
   gPad->SetGrid();
   const int col[NP] = {kGray + 2, kBlue + 2, kRed + 1};
   auto *lg = new TLegend(0.12, 0.72, 0.45, 0.89);
   double ymax = 0;
   for (int p = 0; p < NP; ++p) ymax = std::max(ymax, h[p]->GetMaximum());
   for (int p = 0; p < NP; ++p) {
      TF1 g(Form("g%d", p), "gaus", -1.5, 1.0);
      h[p]->Fit(&g, "QNR");
      printf("  %-30s bound peak %+.3f, sigma %.3f MeV\n", nm[p], g.GetParameter(1), g.GetParameter(2));
      h[p]->SetLineColor(col[p]); h[p]->SetLineWidth(2); h[p]->SetMaximum(1.15 * ymax);
      h[p]->Draw(p ? "hist same" : "hist");
      lg->AddEntry(h[p], nm[p], "l");
   }
   lg->Draw();
   c->SaveAs(png);
}
