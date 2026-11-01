/// @file pid_theta_anatomy_dp.C
/// @brief (theta variant of pid_dedx_anatomy_dp.C) 3D polar-angle estimators, same tracks: tz (production z(phi)
///        line, hit refit off), tc<k> CHORD angle: linear fit of the xy distance from the vertex-end hit vs z over
///        the first f = 25/50/100 % of the window hits (Spyral-like), and rh (hit-refit circle radius, mm).
///        ta<k>: z vs xy ARC LENGTH round the hit-refit circle (slope = cot theta), first 25/50 % of the window hits.
///        tp<k>: MEDIAN of local angles -- 10 mm backbone segments, each segment's 3D principal axis gives
///        theta_i = acos|u_z| (the XZ/YZ slopes of a helix rotate, |d(xy)/dz| = tan theta does not); segments of
///        10 / 25 mm over the whole window. Simulation: simFile + pdg add tst = sin(theta_true) (one-track events only).
///        Kept from the dE/dx anatomy: why is Spyral's dE/dx cleaner than AtHelix3D's for fast tracks? Same tracks (pad-gain corrected), the
///        3D vertex-end segment of firstArcPID rebuilt piece by piece (2026-10-31):
///          d0  production: inlier centroids on the small pads, first inlier -> leaving rho <= 152, / helix path
///          d1  same segment and path, charge of ALL backbone hits in it (no centroid inlier rejection)
///          d2  same segment, production charge / 3D polyline length of the centroids (no helix model)
///          d3  ALL backbone hits / polyline path (d1 charge, d2 path)
///          d4  Spyral-like start: from the vertex-end centroid itself (not the first inlier), all hits, polyline
///          d5  d4 plus the hits blob removal dropped, assigned by z_digi inside the segment's z range
///        Branches as pid_padgain_test_dp.C (sqa/br/sqp, s1/b1, v3 b3 sq3h th3 sp3) plus d0..d5.
///   root -l -b -q 'pid_theta_anatomy_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2_v1.txt","out.root")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
void pid_theta_anatomy_dp(TString in, TString padMap, TString out, Long64_t maxEv = 1000000000, TString simFile = "",
                          int pdg = 0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   TFile *fs = simFile.Length() ? TFile::Open(simFile) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   TClonesArray *mc = nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   std::map<int, double> gain;
   {
      std::ifstream fm(padMap.Length() ? padMap.Data() : "/dev/null");
      std::string line;
      while (std::getline(fm, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         int p; double g; long n;
         if (ss >> p >> g >> n) gain[p] = g;
      }
   }
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   TFile fo(out, "RECREATE");
   TTree o("w", "3D dE/dx anatomy");
   int nh, va, vp, v3, ok1;
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s1, b1, d[6], tz, tch[3], rh, rc, ta[2], tp[2], tst, tr2[2], tq[2];
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s1", &s1); o.Branch("b1", &b1); o.Branch("ok1", &ok1);
   for (int k = 0; k < 6; ++k) o.Branch(Form("d%d", k), &d[k]);
   o.Branch("tz", &tz); o.Branch("rh", &rh); o.Branch("rc", &rc);
   for (int k = 0; k < 3; ++k) o.Branch(Form("tc%d", k), &tch[k]);
   for (int k = 0; k < 2; ++k) { o.Branch(Form("ta%d", k), &ta[k]); o.Branch(Form("tp%d", k), &tp[k]); }
   o.Branch("tst", &tst);
   o.Branch("tr0", &tr2[0]); o.Branch("tr1", &tr2[1]);
   o.Branch("tq0", &tq[0]); o.Branch("tq1", &tq[1]);
   AtTools::Helix3D::EndOpt eoH;
   eoH.hitRefit = true;
   namespace HX = AtTools::Helix3D;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      tst = -1;
      if (ts) {
         int nl = 0;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *mt = (AtMCTrack *)mc->At(j);
            if (!mt || mt->GetPdgCode() != pdg) continue;
            double px = mt->GetPx(), py = mt->GetPy(), pz = mt->GetPz();
            tst = std::hypot(px, py) / std::sqrt(px * px + py * py + pz * pz);
            break;
         }
         if (tst < 0) continue;
      }
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         nh = tr.GetHitArray().size();
         if (nh < 30) continue;
         AtTrack tc(tr);
         for (auto &h : tc.GetHitArray()) {
            auto it = gain.find(h->GetPadNum());
            if (it != gain.end() && it->second > 0) h->SetCharge(h->GetCharge() / it->second);
         }
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
         auto sp = spy.Estimate(tc);
         va = sr.valid; vp = sp.valid;
         br = sr.valid ? sr.brho : -1; sqa = sr.valid ? sr.sqrtdEdx : -1; sqp = sp.valid ? sp.sqrtdEdx : -1;
         pol = sr.valid ? sr.polar * TMath::RadToDeg() : -999;
         auto r0 = h3.Estimate(tr, sr.valid ? sr.direction : -1);
         auto r1 = h3.Estimate(tc, sr.valid ? sr.direction : -1);
         v3 = r0.valid; ok1 = r1.valid;
         b3 = r0.valid ? r0.brho : -1; sq3h = r0.valid ? r0.sqrtdEdx : -1;
         th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         s1 = r1.valid ? r1.sqrtdEdx : -1; b1 = r1.valid ? r1.brho : -1;
         for (auto &x : d) x = -1;
         tz = rh = rc = -1; for (auto &x : tch) x = -1;
         ta[0] = ta[1] = tp[0] = tp[1] = -1;
         tr2[0] = tr2[1] = tq[0] = tq[1] = -1;
         // rebuild the core on the pad-corrected hits (production settings)
         std::vector<HX::P> h;
         for (auto &x : tc.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5);
         if (r.ok) {
            const double zA = h[r.order.front()].z, zB = h[r.order.back()].z;
            int spyEnd = -1;
            if (sr.valid && (sr.direction == 0 || sr.direction == 1)) spyEnd = (sr.direction == 0) == (zA < zB) ? 0 : 1;
            auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
            const int end = V.end;
            const auto &E = r.e[end];
            const int m = r.cen.size();
            if (E.ok && E.brho > 0) {
               tz = E.theta; rc = E.R;
               auto rH = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eoH);
               if (rH.ok && rH.e[end].ok) rh = rH.e[end].R;
               // window hits from the vertex end, lab z
               const double lim = (end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[m - E.nUsed]) + 2.0;
               std::vector<std::pair<double, int>> wv;
               for (size_t i = 0; i < r.order.size(); ++i) {
                  double sf = end == 0 ? r.s[i] : r.length - r.s[i];
                  if (sf <= lim) wv.push_back({sf, r.order[i]});
               }
               std::sort(wv.begin(), wv.end());
               if (wv.size() >= 12) {
                  const auto &h0 = h[wv[0].second];
                  const double fr[3] = {0.25, 0.5, 1.0};
                  for (int k = 0; k < 3; ++k) {
                     int n = std::max(6, (int)(fr[k] * wv.size()));
                     n = std::min<int>(n, wv.size());
                     double sx = 0, sy = 0, sxx = 0, sxy = 0;
                     for (int i = 0; i < n; ++i) {
                        const auto &hp = h[wv[i].second];
                        double zl = 1000.0 - hp.z, rr = std::hypot(hp.x - h0.x, hp.y - h0.y);
                        sx += zl; sy += rr; sxx += zl * zl; sxy += zl * rr;
                     }
                     double den = n * sxx - sx * sx;
                     if (std::abs(den) < 1e-9) continue;
                     double slope = (n * sxy - sx * sy) / den; // d(chord)/dz = tan(theta) (sign: travel in -z -> negative)
                     double th = std::atan(std::abs(slope));
                     tch[k] = slope >= 0 ? th : M_PI - th;
                  }
               }
               // z vs xy arc length round the hit-refit circle, near the vertex
               if (wv.size() >= 12 && rh > 0) {
                  const auto &EH = rH.e[end];
                  const double fr[2] = {0.25, 0.5};
                  for (int k = 0; k < 2; ++k) {
                     int n = std::min<int>(wv.size(), std::max(6, (int)(fr[k] * wv.size())));
                     double prev = std::atan2(h[wv[0].second].y - EH.yc, h[wv[0].second].x - EH.xc), ph = 0;
                     double sx = 0, sy = 0, sxx = 0, sxy = 0;
                     for (int i = 0; i < n; ++i) {
                        const auto &hp = h[wv[i].second];
                        double a = std::atan2(hp.y - EH.yc, hp.x - EH.xc), dph = a - prev;
                        while (dph > M_PI) dph -= 2 * M_PI;
                        while (dph < -M_PI) dph += 2 * M_PI;
                        if (i) ph += dph;
                        prev = a;
                        double sl = EH.R * std::abs(ph), zl = 1000.0 - hp.z;
                        sx += sl; sy += zl; sxx += sl * sl; sxy += sl * zl;
                     }
                     double den = n * sxx - sx * sx;
                     if (std::abs(den) < 1e-9) continue;
                     double slope = (n * sxy - sx * sy) / den; // dz/ds_xy = cot(theta)
                     ta[k] = std::atan2(1.0, slope);
                  }
               }
               // z vs xy arc length ALONG THE TRACK (centroid polyline), first 50 % of the window
               if (E.nUsed >= 6) {
                  auto atq = [&](int jj) { return end == 0 ? jj : m - 1 - jj; };
                  const int nw = std::min(E.nUsed, m);
                  std::vector<double> cs, cl, cz; // backbone s from this end, cumulative xy length, lab z
                  double L = 0;
                  for (int jj = 0; jj < nw; ++jj) {
                     const auto &p = r.cen[atq(jj)];
                     if (jj) { const auto &u = r.cen[atq(jj - 1)]; L += std::hypot(p.x - u.x, p.y - u.y); }
                     cs.push_back(end == 0 ? r.sc[atq(jj)] : r.length - r.sc[atq(jj)]); cl.push_back(L); cz.push_back(p.z);
                  }
                  auto lsq = [&](const std::vector<double> &X, const std::vector<double> &Y, double &th) {
                     const int n = X.size();
                     if (n < 5) return;
                     double sx = 0, sy = 0, sxx = 0, sxy = 0;
                     for (int i = 0; i < n; ++i) { sx += X[i]; sy += Y[i]; sxx += X[i] * X[i]; sxy += X[i] * Y[i]; }
                     double den = n * sxx - sx * sx;
                     if (std::abs(den) > 1e-9) th = std::atan2(1.0, (n * sxy - sx * sy) / den);
                  };
                  const int nh50 = std::max(5, nw / 2);
                  lsq(std::vector<double>(cl.begin(), cl.begin() + nh50), std::vector<double>(cz.begin(), cz.begin() + nh50), tq[0]);
                  // hits: xy length interpolated at their backbone coordinate, those within the first half of the window
                  const double sCut = cs[nh50 - 1];
                  std::vector<double> X, Y;
                  for (auto &pr : wv) {
                     if (pr.first > sCut) break;
                     double sf = pr.first, lx;
                     if (sf <= cs.front()) lx = cl.front() - (cs.front() - sf);
                     else {
                        size_t k = std::upper_bound(cs.begin(), cs.end(), sf) - cs.begin();
                        if (k >= cs.size()) k = cs.size() - 1;
                        double f = (sf - cs[k - 1]) / std::max(1e-9, cs[k] - cs[k - 1]);
                        lx = cl[k - 1] + f * (cl[k] - cl[k - 1]);
                     }
                     X.push_back(lx); Y.push_back(1000.0 - h[pr.second].z);
                  }
                  lsq(X, Y, tq[1]);
               }
               // ROBUST z vs arc length (first 50 %): iterative 2.5 x MAD clipping; tr1 also drops hits with rho < 30 mm
               if (wv.size() >= 12 && rh > 0) {
                  const auto &EH = rH.e[end];
                  int n = std::min<int>(wv.size(), std::max(6, (int)(0.5 * wv.size())));
                  std::vector<double> sA, zA, rhoA;
                  double prev = std::atan2(h[wv[0].second].y - EH.yc, h[wv[0].second].x - EH.xc), ph = 0;
                  for (int i = 0; i < n; ++i) {
                     const auto &hp = h[wv[i].second];
                     double a = std::atan2(hp.y - EH.yc, hp.x - EH.xc), dph = a - prev;
                     while (dph > M_PI) dph -= 2 * M_PI;
                     while (dph < -M_PI) dph += 2 * M_PI;
                     if (i) ph += dph;
                     prev = a;
                     sA.push_back(EH.R * std::abs(ph)); zA.push_back(1000.0 - hp.z); rhoA.push_back(std::hypot(hp.x, hp.y));
                  }
                  for (int v = 0; v < 2; ++v) {
                     std::vector<char> use(sA.size(), 1);
                     if (v == 1) for (size_t i = 0; i < sA.size(); ++i) if (rhoA[i] < 30) use[i] = 0;
                     double slope = 0, icpt = 0;
                     bool okf = false;
                     for (int it = 0; it < 4; ++it) {
                        double sx = 0, sy = 0, sxx = 0, sxy = 0; int nn = 0;
                        for (size_t i = 0; i < sA.size(); ++i) if (use[i]) { sx += sA[i]; sy += zA[i]; sxx += sA[i] * sA[i]; sxy += sA[i] * zA[i]; ++nn; }
                        double den = nn * sxx - sx * sx;
                        if (nn < 6 || std::abs(den) < 1e-9) { okf = false; break; }
                        slope = (nn * sxy - sx * sy) / den; icpt = (sy - slope * sx) / nn; okf = true;
                        std::vector<double> res;
                        for (size_t i = 0; i < sA.size(); ++i) if (use[i]) res.push_back(std::abs(zA[i] - icpt - slope * sA[i]));
                        std::sort(res.begin(), res.end());
                        double mad = std::max(0.5, res[res.size() / 2]);
                        bool changed = false;
                        for (size_t i = 0; i < sA.size(); ++i)
                           if (use[i] && std::abs(zA[i] - icpt - slope * sA[i]) > 2.5 * 1.4826 * mad) { use[i] = 0; changed = true; }
                        if (!changed) break;
                     }
                     if (okf) tr2[v] = std::atan2(1.0, slope);
                  }
               }
               // median of local 3D principal-axis angles over 10 / 25 mm backbone segments
               if (wv.size() >= 12) {
                  const double segL[2] = {10.0, 25.0};
                  const double sEnd = wv.back().first;
                  for (int k = 0; k < 2; ++k) {
                     std::vector<double> ang;
                     size_t i0 = 0;
                     while (i0 < wv.size() && wv[i0].first <= sEnd) {
                        size_t i1 = i0;
                        while (i1 < wv.size() && wv[i1].first - wv[i0].first < segL[k]) ++i1;
                        if (i1 - i0 >= 4) {
                           double mx = 0, my = 0, mz = 0;
                           for (size_t i = i0; i < i1; ++i) { mx += h[wv[i].second].x; my += h[wv[i].second].y; mz += h[wv[i].second].z; }
                           const double nn = i1 - i0;
                           mx /= nn; my /= nn; mz /= nn;
                           TMatrixDSym C(3);
                           for (size_t i = i0; i < i1; ++i) {
                              double dv[3] = {h[wv[i].second].x - mx, h[wv[i].second].y - my, h[wv[i].second].z - mz};
                              for (int a = 0; a < 3; ++a) for (int bq = 0; bq < 3; ++bq) C(a, bq) += dv[a] * dv[bq];
                           }
                           TMatrixDSymEigen eg(C);
                           auto ev = eg.GetEigenVectors(); // column 0 = largest eigenvalue
                           double uz = std::abs(ev(2, 0)), un = std::sqrt(ev(0, 0) * ev(0, 0) + ev(1, 0) * ev(1, 0) + uz * uz);
                           ang.push_back(std::acos(std::min(1.0, uz / un)));
                        }
                        i0 = i1;
                     }
                     if (ang.size() >= 2) { std::sort(ang.begin(), ang.end()); tp[k] = ang[ang.size() / 2]; }
                  }
               }
               // replicate firstArcPID, remembering the segment's centroid range [jS, jE] (traversal index)
               auto inl = [&](const HX::P &p) { return std::abs(std::hypot(p.x - E.xc, p.y - E.yc) - E.R) < E.resCut; };
               auto at = [&](int jj) { return end == 0 ? jj : m - 1 - jj; };
               int j0 = at(0), jFar = 0;
               for (int jj = 0; jj < m; ++jj) if (inl(r.cen[at(jj)])) { j0 = at(jj); break; }
               double dFar = -1;
               for (int jj = 0; jj < m; ++jj) {
                  const auto &p = r.cen[at(jj)];
                  if (!inl(p)) continue;
                  double dd = std::hypot(p.x - r.cen[j0].x, p.y - r.cen[j0].y);
                  if (dd > dFar) { dFar = dd; jFar = jj; }
               }
               double q = 0, phPrev = 0, phAcc = 0, phFirst = 0, phLast = 0;
               bool started = false;
               int jS = -1, jE = -1;
               for (int jj = 0; jj <= jFar; ++jj) {
                  const auto &p = r.cen[at(jj)];
                  double a = std::atan2(p.y - E.yc, p.x - E.xc);
                  if (jj == 0) phPrev = a;
                  else {
                     double dph = a - phPrev;
                     while (dph > M_PI) dph -= 2 * M_PI;
                     while (dph < -M_PI) dph += 2 * M_PI;
                     phAcc += dph; phPrev = a;
                  }
                  bool inner = std::hypot(p.x, p.y) <= 152.0;
                  if (!inl(p)) continue;
                  if (!inner && started) break;
                  if (inner) {
                     if (!started) { phFirst = phAcc; started = true; jS = jj; }
                     phLast = phAcc; jE = jj; q += p.q;
                  }
               }
               const double path = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi) * std::abs(phLast - phFirst);
               if (started && path >= 5) {
                  d[0] = std::sqrt(q / path);
                  // backbone-coordinate range of the segment (centroid sc are from end A)
                  auto sFromEnd = [&](double sA) { return end == 0 ? sA : r.length - sA; };
                  const double sLo = sFromEnd(r.sc[at(jS)]) - 2.0, sHi = sFromEnd(r.sc[at(jE)]) + 2.0;
                  double qAll = 0, qAll0 = 0;
                  double zMin = 1e9, zMax = -1e9;
                  for (size_t i = 0; i < r.order.size(); ++i) {
                     double sf = sFromEnd(r.s[i]);
                     const auto &hp = h[r.order[i]];
                     if (std::hypot(hp.x, hp.y) > 152.0) continue;
                     if (sf >= sLo && sf <= sHi) qAll += hp.q;
                     if (sf >= -1e9 && sf <= sHi) { qAll0 += hp.q; zMin = std::min(zMin, hp.z); zMax = std::max(zMax, hp.z); }
                  }
                  // 3D polyline length of the centroids over the segment (lab z frame like r.cen)
                  auto poly = [&](int ja, int jb) {
                     double L = 0;
                     for (int jj = ja + 1; jj <= jb; ++jj) {
                        const auto &u = r.cen[at(jj - 1)], &v = r.cen[at(jj)];
                        L += std::sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) + (u.z - v.z) * (u.z - v.z));
                     }
                     return L;
                  };
                  const double pl = poly(jS, jE), pl0 = poly(0, jE);
                  d[1] = std::sqrt(qAll / path);
                  if (pl > 5) { d[2] = std::sqrt(q / pl); d[3] = std::sqrt(qAll / pl); }
                  if (pl0 > 5) {
                     d[4] = std::sqrt(qAll0 / pl0);
                     // hits dropped by blob removal (not in r.order) inside the segment's digi-z range, small pads
                     std::vector<char> inOrd(h.size(), 0);
                     for (int i : r.order) inOrd[i] = 1;
                     double qb = 0;
                     for (size_t i = 0; i < h.size(); ++i)
                        if (!inOrd[i] && std::hypot(h[i].x, h[i].y) <= 152.0 && h[i].z >= zMin && h[i].z <= zMax) qb += h[i].q;
                     d[5] = std::sqrt((qAll0 + qb) / pl0);
                  }
               }
            }
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
