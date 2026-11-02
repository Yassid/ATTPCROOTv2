/// @file classify_short_t_dp.C
/// @brief Automatic classification of the data tracks the circle dE/dx PID calls TRITON-LIKE but whose range is shorter
///        than a triton's (ln R/R_t < -0.3), against two controls: triton-like tracks ON the triton range curve
///        (|ln R/R_t| < 0.15) and deuteron-band tracks (2026-11-01). All four runs, circle Brho bLo-bHi, stopping tracks.
///        Per track (hits ordered from the vertex end along the AtHelix3D backbone):
///          MIDPEAK  the charge profile (10 mm bins, 3-bin smoothed) peaks before 75 % of the track, more than 30 mm
///                   before its end, and the mean charge after peak+20 mm is < 35 % of the peak: the particle stopped
///                   and the candidate continues with something else
///          ZREV     z vs path rises then falls (or the reverse): slopes of the first and last 40 % opposite in sign, each
///                   moving > 8 mm in z -- two tracks joined at a vertex
///          WIDE     > 15 % of the hits lie > 10 mm from every backbone centroid: a branch or a second track in the hits
///          STRAIGHT the vertex-end circle turns < 0.5 rad over the track and the range is < 250 mm: a short, nearly
///                   straight stopper whose circle Brho is unreliable
///          CLEAN    none of the above
///        First matching class wins in that order. Also draws nShow examples (xy, saturated hits red; amplitude vs path).
///   root -l -b -q 'classify_short_t_dp.C("plots/short_t_examples.png")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
struct Feat { int cls{-1}; double fpk{0}, tail{0}, turn{0}, wide{0}, L{0}; int zrev{0}; };
Feat classify(const std::vector<AtTools::Helix3D::P> &h, const AtTools::Helix3D::Result &r, int end, std::vector<double> &S,
              std::vector<int> &O)
{
   Feat F;
   const int n = r.order.size();
   std::vector<std::pair<double, int>> ord;
   for (int i = 0; i < n; ++i) ord.push_back({end == 0 ? r.s[i] : r.length - r.s[i], r.order[i]});
   std::sort(ord.begin(), ord.end());
   S.clear(); O.clear();
   for (auto &o : ord) { S.push_back(o.first); O.push_back(o.second); }
   F.L = r.length;
   // charge profile
   const int nb = std::max(1, (int)std::ceil(F.L / 10.0));
   std::vector<double> prof(nb, 0), sm(nb, 0);
   for (int i = 0; i < n; ++i) prof[std::min(nb - 1, (int)(S[i] / 10.0))] += h[O[i]].q;
   for (int b = 0; b < nb; ++b) { double s = 0; int c = 0; for (int d = -1; d <= 1; ++d) if (b + d >= 0 && b + d < nb) { s += prof[b + d]; ++c; } sm[b] = s / c; }
   int bp = std::max_element(sm.begin(), sm.end()) - sm.begin();
   const double sPk = (bp + 0.5) * 10;
   F.fpk = sPk / std::max(1.0, F.L);
   double ta = 0; int tc = 0;
   for (int b = bp + 2; b < nb; ++b) { ta += sm[b]; ++tc; }
   F.tail = tc ? ta / tc / std::max(1e-9, sm[bp]) : 1;
   const bool midpeak = F.fpk < 0.75 && F.L - sPk > 30 && F.tail < 0.35;
   // z reversal
   auto slope = [&](int i0, int i1, double &dz) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0; int m = 0;
      for (int i = i0; i < i1; ++i) { sx += S[i]; sy += h[O[i]].z; sxx += S[i] * S[i]; sxy += S[i] * h[O[i]].z; ++m; }
      double den = m * sxx - sx * sx;
      double b = std::abs(den) > 1e-9 ? (m * sxy - sx * sy) / den : 0;
      dz = b * (S[i1 - 1] - S[i0]);
      return b;
   };
   double dz1, dz2;
   double b1 = slope(0, std::max(3, (int)(0.4 * n)), dz1), b2 = slope(std::min(n - 3, (int)(0.6 * n)), n, dz2);
   F.zrev = (b1 * b2 < 0 && std::abs(dz1) > 8 && std::abs(dz2) > 8);
   // width: hits far from every backbone centroid
   int far = 0;
   for (int i = 0; i < n; ++i) {
      double dmin = 1e9;
      // AtHelix3D's centroids live in its internal frame z' = 1000 - z (zSign -1, zOff 1000): compare like with like
      const double zi = 1000.0 - h[i].z;
      for (auto &c : r.cen) dmin = std::min(dmin, std::sqrt((h[i].x - c.x) * (h[i].x - c.x) + (h[i].y - c.y) * (h[i].y - c.y) + (zi - c.z) * (zi - c.z)));
      far += dmin > 10;
   }
   F.wide = double(far) / n;
   // turning of the vertex-end circle along the centroids
   const auto &E = r.e[end];
   double prev = 0, acc = 0;
   const int m = r.cen.size();
   for (int j = 0; j < m; ++j) {
      const auto &c = r.cen[end == 0 ? j : m - 1 - j];
      double a = std::atan2(c.y - E.yc, c.x - E.xc);
      if (j) { double d = a - prev; while (d > M_PI) d -= 2 * M_PI; while (d < -M_PI) d += 2 * M_PI; acc += d; }
      prev = a;
   }
   F.turn = std::abs(acc);
   F.cls = midpeak ? 0 : F.zrev ? 1 : F.wide > 0.15 ? 2 : (F.turn < 0.5 && F.L < 250) ? 3 : 4;
   return F;
}
void classify_short_t_dp(TString png, double bLo = 0.15, double bHi = 0.35, int nShow = 12,
                         TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root", TString rr = "/mnt/f/a1975/pid_range/R_all.root")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   namespace HX = AtTools::Helix3D;
   std::vector<double> bx, rt;
   {
      TChain t("w"), tr("r");
      for (const char *s : {"gs", "ex2"}) { t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); tr.Add(Form("/mnt/f/a1975/pid_range/r_%s.root", s)); }
      t.AddFriend(&tr, "r");
      for (double b = 0.10; b < 0.48; b += 0.02) {
         long n = t.Draw("r.len", Form("va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<955&&bA>%g&&bA<%g", b, b + 0.02), "goff");
         if (n < 8) continue;
         std::vector<double> v(t.GetV1(), t.GetV1() + n);
         std::sort(v.begin(), v.end());
         bx.push_back(b + 0.01); rt.push_back(v[n / 2]);
      }
   }
   auto rtAt = [&](double b) { for (size_t i = 0; i + 1 < bx.size(); ++i) if (b >= bx[i] && b < bx[i + 1]) return rt[i] + (rt[i + 1] - rt[i]) * (b - bx[i]) / (bx[i + 1] - bx[i]); return -1.; };
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("r=r", rr);
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<955&&r.len>20";
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   // deuteron mode per 0.05 slice
   std::vector<double> dmv;
   for (double b = 0.15; b < 0.351; b += 0.05) {
      TH1D hh("hh", "", 600, -4, 2);
      w->Draw(Form("log(s1_g/%g*bA)>>hh", md), sel + Form("&&bA>%g&&bA<%g", b, b + 0.05), "goff");
      hh.Smooth(2); dmv.push_back(hh.GetBinCenter(hh.GetMaximumBin()));
   }
   w->SetEstimate(w->GetEntries() + 1);
   long n = w->Draw(Form("Entry$:bA:r.len:log(s1_g/%g*bA)", md), sel + Form("&&bA>%g&&bA<%g", bLo, bHi), "goff");
   // group: 0 short t-like, 1 on-curve t-like, 2 deuteron band (every 20th, a control sample)
   std::map<Long64_t, int> grp;
   long nd = 0;
   for (long i = 0; i < n; ++i) {
      double b = w->GetV2()[i];
      int k = std::min<int>(dmv.size() - 1, std::max(0, (int)((b - 0.15) / 0.05)));
      double lnm = w->GetV4()[i] - dmv[k], r0 = rtAt(b);
      if (r0 <= 0) continue;
      double lrt = std::log(w->GetV3()[i] / r0);
      Long64_t row = (Long64_t)w->GetV1()[i];
      if (lnm > 0.25 && lnm < 0.55) { if (lrt < -0.3) grp[row] = 0; else if (std::abs(lrt) < 0.15) grp[row] = 1; }
      else if (std::abs(lnm) < 0.12 && (nd++ % 20) == 0) grp[row] = 2;
   }
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   const char *cn[5] = {"MIDPEAK", "ZREV", "WIDE", "STRAIGHT", "CLEAN"};
   long cnt[3][5] = {{0}};
   struct Ex { TString run; Long64_t ev; int itr, cls; double L, lrt, bA; };
   std::vector<Ex> ex;
   std::map<Long64_t, double> rowB, rowL;
   for (long i = 0; i < n; ++i) { rowB[(Long64_t)w->GetV1()[i]] = w->GetV2()[i]; rowL[(Long64_t)w->GetV1()[i]] = w->GetV3()[i]; }
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   HX::EndOpt eo;
   eo.hitRefit = true; eo.hitRefitZ = true;
   Long64_t off = 0;
   for (int ir = 0; ir < 4; ++ir) {
      TFile fr(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", runs[ir]));
      auto *t = (TTree *)fr.Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      Long64_t row = off;
      for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
         if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
         int itr = 0;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            if (tr.GetHitArray().size() >= 30) {
               auto it = grp.find(row);
               if (it != grp.end()) {
                  std::vector<HX::P> h;
                  for (auto &q : tr.GetHitArray()) { auto p = q->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), q->GetCharge()}); }
                  auto sr = spy.Estimate(tr);
                  auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
                  if (r.ok) {
                     int spyEnd = (sr.valid && (sr.direction == 0 || sr.direction == 1)) ? ((sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1) : -1;
                     auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
                     std::vector<double> S; std::vector<int> O;
                     auto F = classify(h, r, V.end, S, O);
                     cnt[it->second][F.cls]++;
                     if (it->second == 0) ex.push_back({runs[ir], ev, itr, F.cls, rowL[row], std::log(rowL[row] / rtAt(rowB[row])), rowB[row]});
                  }
               }
               ++row;
            }
            ++itr;
         }
      }
      off = row;
   }
   const char *gn[3] = {"SHORT triton-like", "triton-like ON t curve", "deuteron band (1/20)"};
   printf("circle Brho %.2f-%.2f, stopping tracks, all four runs\n%-24s %6s", bLo, bHi, "group", "N");
   for (int c = 0; c < 5; ++c) printf(" %9s", cn[c]);
   printf("\n");
   for (int g = 0; g < 3; ++g) {
      long tot = 0; for (int c = 0; c < 5; ++c) tot += cnt[g][c];
      printf("%-24s %6ld", gn[g], tot);
      for (int c = 0; c < 5; ++c) printf("  %6.1f%%", 100. * cnt[g][c] / std::max(tot, 1L));
      printf("\n");
   }
   // examples: spread over the classes
   std::vector<Ex> show;
   for (int c = 0; c < 5 && (int)show.size() < nShow; ++c) {
      int k = 0;
      for (auto &e : ex) if (e.cls == c && k < std::max(2, nShow / 5 + 1) && (int)show.size() < nShow) { show.push_back(e); ++k; }
   }
   const int ncol = 4, nrow = (show.size() + 1) / 2;
   auto *cv = new TCanvas("cv", "cv", 1800, 380 * nrow);
   cv->Divide(ncol, nrow);
   for (size_t k = 0; k < show.size(); ++k) {
      auto &e = show[k];
      TFile fr(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", e.run.Data()));
      auto *t = (TTree *)fr.Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      t->GetEntry(e.ev);
      auto &tr = ((AtPatternEvent *)pa->At(0))->GetTrackCand()[e.itr];
      std::vector<HX::P> h;
      for (auto &q : tr.GetHitArray()) { auto p = q->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), q->GetCharge()}); }
      auto sr = spy.Estimate(tr);
      auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
      int spyEnd = (sr.valid && (sr.direction == 0 || sr.direction == 1)) ? ((sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1) : -1;
      auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
      std::vector<double> S; std::vector<int> O;
      auto F = classify(h, r, V.end, S, O);
      std::vector<double> X, Y, Q, xs, ys;
      for (size_t i = 0; i < O.size(); ++i) { X.push_back(h[O[i]].x); Y.push_back(h[O[i]].y); Q.push_back(h[O[i]].q); if (h[O[i]].q > 3500) { xs.push_back(X.back()); ys.push_back(Y.back()); } }
      cv->cd(2 * k + 1);
      auto *g = new TGraph(X.size(), X.data(), Y.data());
      g->SetTitle(Form("%s: run %s ev %lld trk %d, B#rho %.2f;x [mm];y [mm]", cn[F.cls], e.run.Data(), e.ev, e.itr, e.bA));
      g->SetMarkerStyle(20); g->SetMarkerSize(0.4); g->SetMarkerColor(kGray + 2); g->Draw("AP");
      if (!xs.empty()) { auto *gs = new TGraph(xs.size(), xs.data(), ys.data()); gs->SetMarkerStyle(20); gs->SetMarkerSize(0.5); gs->SetMarkerColor(kRed); gs->Draw("P"); }
      auto *mv = new TMarker(X[0], Y[0], 29); mv->SetMarkerColor(kGreen + 2); mv->SetMarkerSize(2); mv->Draw();
      cv->cd(2 * k + 2);
      auto *gq = new TGraph(S.size(), S.data(), Q.data());
      gq->SetTitle(Form("range %.0f (R/R_{t} %.2f) peak at %.0f%%, tail %.2f, turn %.2f rad, wide %.2f;path from vertex end [mm];amplitude [ADC]",
                        e.L, std::exp(e.lrt), 100 * F.fpk, F.tail, F.turn, F.wide));
      gq->SetMarkerStyle(20); gq->SetMarkerSize(0.4); gq->Draw("AP");
      auto *l = new TLine(0, 3500, S.back(), 3500); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw();
   }
   cv->SaveAs(png);
}
