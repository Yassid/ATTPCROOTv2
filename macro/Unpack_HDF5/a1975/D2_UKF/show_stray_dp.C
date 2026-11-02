/// @file show_stray_dp.C
/// @brief Event displays of stray hits (AtHelix3D blob removal), 2026-11-02: hits kept grey, dropped red, hits the large-pad guard keeps
///        orange (dropped by the removal but beyond 152 mm). Four cases taken from the order_fit_dp.C results: run 0020 event 9789; a data
///        deuteron-band track from run 0026 that production gets wrong and the guard recovers; a fast simulated deuteron likewise; a slow simulated
///        deuteron whose real outer track the unguarded removal cuts (mode 8 wrong, mode 9 right). xy view, and z_digi vs radius from the axis.
///   root -l -b -q 'show_stray_dp.C("plots/stray_examples.png")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
namespace strayshow {
const double kMd = 1875.613;
double br(double ke) { return std::sqrt(ke * ke + 2 * ke * kMd) / 299.792458; }
double exEl(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ed = ke + kMd, pd = std::sqrt(Ed * Ed - kMd * kMd);
   const double E4 = Eb + kMd - Ed, px = -pd * std::sin(th), pz = pb - pd * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - mC;
}
/// the n-th >=30-hit track of the file (sim: only events with exactly one such track and a truth particle)
bool findRow(TTree *t, TClonesArray *&pa, TTree *ts, TClonesArray *&mc, int pdg, Long64_t row, AtTrack &out, Long64_t &evOut)
{
   Long64_t count = 0;
   for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto &c = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      if (ts) {
         int nl = 0;
         for (auto &tr : c) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         bool has = false;
         for (int j = 0; j < mc->GetEntriesFast(); j++) { auto *m = (AtMCTrack *)mc->At(j); if (m && m->GetPdgCode() == pdg) has = true; }
         if (!has) continue;
      }
      for (auto &tr : c) {
         if (tr.GetHitArray().size() < 30) continue;
         if (count == row) { out = tr; evOut = ev; return true; }
         ++count;
      }
   }
   return false;
}
void draw(AtTrack &tr, TCanvas *c, int rowPad, TString title)
{
   namespace HX = AtTools::Helix3D;
   std::vector<HX::P> h;
   for (auto &x : tr.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
   HX::EndOpt eo;
   auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
   std::set<int> kept(r.order.begin(), r.order.end());
   std::vector<double> X[3], Y[3], Z[3], R[3];
   for (size_t i = 0; i < h.size(); ++i) {
      const double rho = std::hypot(h[i].x, h[i].y);
      int k = kept.count(i) ? 0 : (rho > 152 ? 2 : 1);
      X[k].push_back(h[i].x); Y[k].push_back(h[i].y); Z[k].push_back(h[i].z); R[k].push_back(rho);
   }
   const int col[3] = {kGray + 1, kRed, kOrange + 1};
   const double sz[3] = {0.45, 0.9, 0.9};
   for (int p = 0; p < 2; ++p) {
      c->cd(2 * rowPad + 1 + p);
      gPad->SetGrid();
      std::vector<double> &a0 = p == 0 ? X[0] : R[0], &b0 = p == 0 ? Y[0] : Z[0];
      auto *g0 = new TGraph(a0.size(), a0.data(), b0.data());
      g0->SetTitle(p == 0 ? Form("%s;x [mm];y [mm]", title.Data()) : Form("%zu kept, %zu dropped (red), %zu kept by the guard (orange);radius from axis [mm];z_{digi} [mm]", X[0].size(), X[1].size(), X[2].size()));
      g0->SetMarkerStyle(20); g0->SetMarkerSize(sz[0]); g0->SetMarkerColor(col[0]);
      // frame on all hits
      double xl = 1e9, xh = -1e9, yl = 1e9, yh = -1e9;
      for (int k = 0; k < 3; ++k) {
         std::vector<double> &a = p == 0 ? X[k] : R[k], &b = p == 0 ? Y[k] : Z[k];
         for (size_t i = 0; i < a.size(); ++i) { xl = std::min(xl, a[i]); xh = std::max(xh, a[i]); yl = std::min(yl, b[i]); yh = std::max(yh, b[i]); }
      }
      g0->GetXaxis()->SetLimits(xl - 0.05 * (xh - xl + 1), xh + 0.05 * (xh - xl + 1));
      g0->SetMinimum(yl - 0.05 * (yh - yl + 1)); g0->SetMaximum(yh + 0.05 * (yh - yl + 1));
      g0->Draw("AP");
      for (int k = 1; k < 3; ++k) {
         std::vector<double> &a = p == 0 ? X[k] : R[k], &b = p == 0 ? Y[k] : Z[k];
         if (a.empty()) continue;
         auto *g = new TGraph(a.size(), a.data(), b.data());
         g->SetMarkerStyle(20); g->SetMarkerSize(sz[k]); g->SetMarkerColor(col[k]);
         g->Draw("P");
      }
   }
}
}
void show_stray_dp(TString png, TString O = "/mnt/f/a1975/order_test/")
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kWarning; gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
   auto *c = new TCanvas("c", "c", 1400, 1760);
   c->Divide(2, 4);
   // --- 1. run 0020 event 9789
   {
      TFile f("/mnt/f/a1975/reco_d2_dv1104/run_0020_multifit_reco.root"); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      t->GetEntry(9789);
      AtTrack tr;
      for (auto &x : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) if (x.GetHitArray().size() >= 30) tr = x;
      strayshow::draw(tr, c, 0, "data run 0020 ev 9789, proton: 0.43 MeV, 83#circ, vertex 83 mm off #rightarrow 3.18 MeV, 73#circ, 8 mm");
   }
   // --- 2. data deuteron: production |Ex| > 2, guard |Ex| < 0.5
   {
      TChain a("o"), g("o");
      a.Add(O + "data0.root"); a.Add(O + "data1.root"); g.Add(O + "g_data0.root"); g.Add(O + "g_data1.root");
      double ke0, th0, ke9, th9, c20, nd0, c29, nd9; int ok0, ok9, nDrop, nh;
      a.SetBranchAddress("ke0", &ke0); a.SetBranchAddress("th0", &th0); a.SetBranchAddress("ok0", &ok0); a.SetBranchAddress("c20", &c20); a.SetBranchAddress("nd0", &nd0); a.SetBranchAddress("nh", &nh);
      g.SetBranchAddress("ke9", &ke9); g.SetBranchAddress("th9", &th9); g.SetBranchAddress("ok9", &ok9); g.SetBranchAddress("nDrop", &nDrop); g.SetBranchAddress("c29", &c29); g.SetBranchAddress("nd9", &nd9);
      Long64_t pick = -1; double e0 = 0, e9 = 0, kA = 0, kB = 0;
      for (Long64_t i = 0; i < a.GetEntries() && pick < 0; ++i) {
         a.GetEntry(i); g.GetEntry(i);
         if (!ok0 || !ok9 || nDrop < 5 || nDrop > 40 || nh > 250) continue;
         e0 = strayshow::exEl(ke0, th0, 180); e9 = strayshow::exEl(ke9, th9, 180);
         if (std::abs(e0) > 3 && std::abs(e9) < 0.3 && ke9 < 2) { pick = i; kA = ke0; kB = ke9; }
      }
      if (pick >= 0) {
         TFile f("/mnt/f/a1975/reco_d2_dv1104/run_0026_multifit_reco.root"); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *pa = nullptr, *mc = nullptr;
         t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
         AtTrack tr; Long64_t ev;
         strayshow::findRow(t, pa, nullptr, mc, 0, pick, tr, ev);
         strayshow::draw(tr, c, 1, Form("data run 0026 ev %lld, as deuteron: Ex %+.1f MeV (%.2f MeV) #rightarrow %+.2f MeV (%.2f MeV)", ev, e0, kA, e9, kB));
      }
   }
   // --- 3. fast sim deuteron, production wrong; 4. slow sim deuteron, unguarded removal wrong
   for (int cs = 0; cs < 2; ++cs) {
      TFile fa(O + "d5001.root"), fe(O + "e_d5001.root"), fg(O + "g_d5001.root");
      auto *a = (TTree *)fa.Get("o"); auto *e = (TTree *)fe.Get("o"); auto *g = (TTree *)fg.Get("o");
      double tbr, ke0, ke8, ke9; int ok0, ok8, ok9, nDrop;
      a->SetBranchAddress("tbr", &tbr); a->SetBranchAddress("ke0", &ke0); a->SetBranchAddress("ok0", &ok0);
      e->SetBranchAddress("ke8", &ke8); e->SetBranchAddress("ok8", &ok8); e->SetBranchAddress("nDrop", &nDrop);
      g->SetBranchAddress("ke9", &ke9); g->SetBranchAddress("ok9", &ok9);
      Long64_t pick = -1; double bA = 0, bB = 0, bT = 0;
      for (Long64_t i = 0; i < a->GetEntries() && pick < 0; ++i) {
         a->GetEntry(i); e->GetEntry(i); g->GetEntry(i);
         if (nDrop <= 0 || !ok9 || !ok0) continue;
         const double l0 = std::log(strayshow::br(ke0) / tbr), l9 = std::log(strayshow::br(ke9) / tbr), l8 = ok8 ? std::log(strayshow::br(ke8) / tbr) : 9;
         if (cs == 0 && tbr > 0.6 && tbr < 0.8 && std::abs(l0) > 0.3 && std::abs(l9) < 0.03) { pick = i; bA = strayshow::br(ke0); bB = strayshow::br(ke9); bT = tbr; }
         if (cs == 1 && tbr > 0.28 && tbr < 0.36 && std::abs(l0) < 0.03 && std::abs(l8) > 0.3 && std::abs(l9) < 0.03) { pick = i; bA = ok8 ? strayshow::br(ke8) : -1; bB = strayshow::br(ke9); bT = tbr; }
      }
      if (pick < 0) { printf("case %d: no track found\n", cs); continue; }
      TFile f("/mnt/f/a1975_C16_dd_sim/el_s5001_reco.root"), fs("/mnt/f/a1975_C16_dd_sim/el_s5001_sim.root");
      auto *t = (TTree *)f.Get("cbmsim"); auto *ts = (TTree *)fs.Get("cbmsim"); TClonesArray *pa = nullptr, *mc = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc);
      AtTrack tr; Long64_t ev;
      strayshow::findRow(t, pa, ts, mc, 1000010020, pick, tr, ev);
      strayshow::draw(tr, c, 2 + cs, cs == 0 ? Form("sim d ev %lld, true B#rho %.2f: production %.2f #rightarrow removal %.2f", ev, bT, bA, bB)
                                              : Form("sim d ev %lld, true B#rho %.2f: unguarded removal %.2f, with guard %.2f", ev, bT, bA, bB));
   }
   c->SaveAs(png);
}
