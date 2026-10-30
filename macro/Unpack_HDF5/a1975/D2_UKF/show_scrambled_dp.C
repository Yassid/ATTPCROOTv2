/// @file show_scrambled_dp.C
/// @brief Picture what "scrambled" means for a1975 (d,p) data tracks: the same hits joined in the
///        three orders the reconstruction uses, one PNG per track. No GUI.
///
///  pad 1  HITS in the stored hit-array order       -- what AtPRA's arc-length theta fit walks
///  pad 2  CLUSTERS sorted by z                     -- what AtGenfitter feeds the Kalman filter
///  pad 3  HITS in 3D minimum-spanning-tree order   -- helix3d.h, uses no z and no input order
///  pad 4  z_lab along each sequence (normalised index); a monotonic curve = an ordered track
/// Segments are coloured by position in the sequence (blue first -> red last), so a scrambled
/// order shows as long criss-crossing chords and colours jumping around the ring.
/// Scramble metric in the title: sequence path length / path length of the 3D-MST order of the SAME
/// points (1.0 = as ordered as the MST walk; the cluster ratio uses the clusters' own MST walk).
///
///   root -l -b -q 'show_scrambled_dp.C("<outdir>")'
#include "../UKF/pid/helix3d.h"

namespace {
double pathLen(const std::vector<std::array<double, 3>> &p)
{
   double s = 0;
   for (size_t i = 1; i < p.size(); ++i)
      s += std::sqrt(std::pow(p[i][0] - p[i - 1][0], 2) + std::pow(p[i][1] - p[i - 1][1], 2) +
                     std::pow(p[i][2] - p[i - 1][2], 2));
   return s;
}
void drawSeq(const std::vector<std::array<double, 3>> &p, bool markers)
{
   const int n = p.size();
   for (int i = 1; i < n; ++i) {
      double f = double(i) / n;
      auto *l = new TLine(p[i - 1][0], p[i - 1][1], p[i][0], p[i][1]);
      l->SetLineColor(TColor::GetColor(int(30 + 220 * f), 40, int(250 - 220 * f)));
      l->SetLineWidth(1);
      l->Draw();
   }
   if (markers) {
      auto *g = new TGraph;
      for (auto &q : p) g->AddPoint(q[0], q[1]);
      g->SetMarkerStyle(20); g->SetMarkerSize(0.35); g->SetMarkerColor(kBlack);
      g->Draw("P");
   }
}
} // namespace

void show_scrambled_dp(TString outDir = "/mnt/c/Users/Yassid/Desktop/scrambled_tracks_2026-10-30",
                       TString reco = "/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root", Long64_t maxEv = 6000,
                       int minHits = 150)
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kError;
   gStyle->SetOptStat(0);
   gSystem->mkdir(outDir, kTRUE);
   TFile fr(reco);
   auto *t = (TTree *)fr.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);

   // pass 1: scramble ratio of the stored hit order for every long track
   struct C { Long64_t ev; int tid; double rHit, rClz; int n; };
   std::vector<C> all;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      t->GetEntry(ev);
      if (!pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &hits = tr.GetHitArray();
         if ((int)hits.size() < minHits) continue;
         std::vector<h3d::P> h;
         std::vector<std::array<double, 3>> st;
         for (auto &x : hits) {
            auto p = x->GetPosition();
            h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()});
            st.push_back({p.X(), p.Y(), 1000 - p.Z()});
         }
         std::vector<int> ord; int nb; double len;
         if (!h3d::backbone(h, ord, nb, len) || len <= 0) continue;
         std::vector<std::array<double, 3>> mst;
         for (int i : ord) mst.push_back(st[i]);
         std::vector<std::array<double, 3>> cz;
         std::vector<h3d::P> hcz;
         if (auto *hc = tr.GetHitClusterArray())
            for (auto &c : *hc) { auto p = c.GetPosition(); cz.push_back({p.X(), p.Y(), 1000 - p.Z()}); hcz.push_back({p.X(), p.Y(), p.Z(), 1}); }
         std::sort(cz.begin(), cz.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
         double rc = -1;
         std::vector<int> oc; int nbc; double lc;
         if (cz.size() > 5 && h3d::backbone(hcz, oc, nbc, lc)) {
            std::vector<std::array<double, 3>> cm;
            for (int i : oc) cm.push_back({hcz[i].x, hcz[i].y, 1000 - hcz[i].z});
            rc = pathLen(cz) / pathLen(cm);
         }
         all.push_back({ev, tr.GetTrackID(), pathLen(st) / pathLen(mst), rc, (int)hits.size()});
      }
   }
   std::vector<double> r;
   for (auto &c : all) r.push_back(c.rHit);
   std::sort(r.begin(), r.end());
   printf("%zu tracks >= %d hits; stored-hit-order path / MST-order path: p10 %.2f median %.2f p90 %.2f\n", all.size(),
          minHits, r[r.size() / 10], r[r.size() / 2], r[9 * r.size() / 10]);

   // pick by QUANTILE of the hit-order ratio: lowest 10 %, around the median, top 10 %
   auto q = [&](double f) { return r[std::min(r.size() - 1, size_t(f * r.size()))]; };
   auto pick = [&](double lo, double hi, int k) {
      std::vector<C> s;
      for (auto &c : all) if (c.rHit >= lo && c.rHit < hi) s.push_back(c);
      std::sort(s.begin(), s.end(), [](auto &a, auto &b) { return a.n < b.n; });
      std::vector<C> o;
      for (int i = 0; i < k && !s.empty(); ++i) o.push_back(s[(i + 0.5) * s.size() / k]);
      return o;
   };
   std::vector<std::pair<std::string, std::vector<C>>> groups = {
      {"A_least_scrambled", pick(0, q(0.10), 2)}, {"B_median", pick(q(0.45), q(0.55), 3)},
      {"C_most_scrambled", pick(q(0.90), 1e9, 3)}};

   FILE *idx = fopen(Form("%s/README.txt", outDir.Data()), "w");
   fprintf(idx, "a1975 16C(d,p) DATA, run_0016 (reco_d2_dv1104, the production the adopted find805 fits use)\n"
                "%zu tracks with >= %d hits; stored-hit-order path / MST-order path: p10 %.2f median %.2f p90 %.2f\n\n"
                "Each PNG = one track, the SAME hits joined in three orders:\n"
                "  1 stored hit-array order (AtPRA arc-length theta fit walks this)\n"
                "  2 clusters sorted by z (AtGenfitter's measurement order)\n"
                "  3 3D minimum-spanning-tree order (helix3d.h; no z, no input order)\n"
                "  4 z_lab along each sequence vs normalised index\n"
                "Colour = position in the sequence, blue first -> red last. Ratio = path / path of the 3D-MST order of the same points (1 = as ordered as the MST walk).\n\n",
           all.size(), minHits, r[r.size() / 10], r[r.size() / 2], r[9 * r.size() / 10]);

   for (auto &[gname, list] : groups) {
      for (auto &cc : list) {
         t->GetEntry(cc.ev);
         AtTrack *me = nullptr;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand())
            if (tr.GetTrackID() == cc.tid) me = &tr;
         if (!me) continue;
         std::vector<h3d::P> h;
         std::vector<std::array<double, 3>> st, cz, mst;
         for (auto &x : me->GetHitArray()) {
            auto p = x->GetPosition();
            h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()});
            st.push_back({p.X(), p.Y(), 1000 - p.Z()});
         }
         for (auto &c : *me->GetHitClusterArray()) { auto p = c.GetPosition(); cz.push_back({p.X(), p.Y(), 1000 - p.Z()}); }
         std::sort(cz.begin(), cz.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
         std::vector<int> ord; int nb; double len;
         h3d::backbone(h, ord, nb, len);
         for (int i : ord) mst.push_back(st[i]);

         double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9, zlo = 1e9, zhi = -1e9;
         for (auto &q : st) {
            xlo = std::min(xlo, q[0]); xhi = std::max(xhi, q[0]); ylo = std::min(ylo, q[1]); yhi = std::max(yhi, q[1]);
            zlo = std::min(zlo, q[2]); zhi = std::max(zhi, q[2]);
         }
         double cx = 0.5 * (xlo + xhi), cy = 0.5 * (ylo + yhi), half = 0.55 * std::max(xhi - xlo, yhi - ylo) + 5;
         auto *c = new TCanvas("c", "c", 1800, 520);
         c->Divide(4, 1);
         const char *ttl[3] = {"stored HIT order (AtPRA)", "CLUSTERS by z (genfit)", "3D MST order (helix3d)"};
         std::vector<std::array<double, 3>> *seq[3] = {&st, &cz, &mst};
         for (int k = 0; k < 3; ++k) {
            c->cd(k + 1);
            gPad->DrawFrame(cx - half, cy - half, cx + half, cy + half,
                            Form("%s   ratio %.2f;x [mm];y [mm]", ttl[k], k == 0 ? cc.rHit : (k == 1 ? cc.rClz : 1.0)));
            drawSeq(*seq[k], k != 1);
            auto *ax = new TMarker(0, 0, 5); ax->SetMarkerSize(1.5); ax->Draw();
         }
         c->cd(4);
         gPad->DrawFrame(0, zlo - 20, 1, zhi + 20, "z_{lab} along the sequence;normalised index;z_{lab} [mm]");
         int col[3] = {kGray + 1, kOrange + 7, kBlue + 1};
         auto *lg = new TLegend(0.12, 0.75, 0.6, 0.89);
         for (int k = 0; k < 3; ++k) {
            auto *g = new TGraph;
            for (size_t i = 0; i < seq[k]->size(); ++i) g->AddPoint(double(i) / seq[k]->size(), (*seq[k])[i][2]);
            g->SetMarkerStyle(20); g->SetMarkerSize(k == 0 ? 0.25 : 0.4); g->SetMarkerColor(col[k]);
            g->Draw("P");
            lg->AddEntry(g, ttl[k], "p");
         }
         lg->Draw();
         TString name = Form("%s_ev%lld_tid%d_n%d.png", gname.c_str(), cc.ev, cc.tid, cc.n);
         c->SaveAs(outDir + "/" + name);
         fprintf(idx, "%-46s hits %4d  hit-order ratio %6.2f  z-sorted-cluster ratio %6.2f\n", name.Data(), cc.n,
                 cc.rHit, cc.rClz);
         delete c;
      }
   }
   fclose(idx);
}
