/// @file pid_dedx_anatomy_dp.C
/// @brief Why is Spyral's dE/dx cleaner than AtHelix3D's for fast tracks? Same tracks (pad-gain corrected), the
///        3D vertex-end segment of firstArcPID rebuilt piece by piece (2026-10-31):
///          d0  production: inlier centroids on the small pads, first inlier -> leaving rho <= 152, / helix path
///          d1  same segment and path, charge of ALL backbone hits in it (no centroid inlier rejection)
///          d2  same segment, production charge / 3D polyline length of the centroids (no helix model)
///          d3  ALL backbone hits / polyline path (d1 charge, d2 path)
///          d4  Spyral-like start: from the vertex-end centroid itself (not the first inlier), all hits, polyline
///          d5  d4 plus the hits blob removal dropped, assigned by z_digi inside the segment's z range
///        Branches as pid_padgain_test_dp.C (sqa/br/sqp, s1/b1, v3 b3 sq3h th3 sp3) plus d0..d5.
///   root -l -b -q 'pid_dedx_anatomy_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2_v1.txt","out.root")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
void pid_dedx_anatomy_dp(TString in, TString padMap, TString out, Long64_t maxEv = 1000000000)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   std::map<int, double> gain;
   {
      std::ifstream fm(padMap.Data());
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
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s1, b1, d[6];
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s1", &s1); o.Branch("b1", &b1); o.Branch("ok1", &ok1);
   for (int k = 0; k < 6; ++k) o.Branch(Form("d%d", k), &d[k]);
   namespace HX = AtTools::Helix3D;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
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
