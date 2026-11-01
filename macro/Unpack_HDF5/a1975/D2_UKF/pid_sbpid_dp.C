/// @file pid_sbpid_dp.C
/// @brief SPECIES-BLIND PID and Brho, all tracks (2026-11-01): pid_bragg_dp.C (arc-length angle, Bragg-curve fit of
///        tracks that stop on the small pads) plus two species-blind equation-of-motion paths (eom_brho_dp.h,
///        velocity-scaled proton CATIMA stopping, Z = 1):
///          through-going: mass from the measured mean dE/dx of the PID segment (solveMass), iterated with the planar
///                         EOM momentum fit (theta = arc-length angle)  -> ef, em  (as pid_dedxeom_dp.C)
///          stopping:      every contained track with a converged Bragg fit: mass from the Bragg kinetic energy and the
///                         momentum, m = (p^2 - KE^2) / 2 KE, iterated with the EOM momentum fit -> ef2, em2
///        Combined: emc / efc = stopping branch if the Bragg rise brr > braggMin or the dE/dx solver failed (no mass or
///        railed at its 8000 MeV limit), else the through-going one; ksel = 1 through, 2 stopping.
///        k = per-run MeV per charge unit from dedxeom_calib_dp.C WITH the EOM momentum (eomPrefix).
///   root -l -b -q 'pid_sbpid_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0026_multifit_reco.root","pid/padgain_d2_v1.txt","out.root",4.97805e-06)'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
#include "eom_brho_dp.h"
#include "bragg_fit_dp.h"
void pid_sbpid_dp(TString in, TString padMap, TString out, double kcal, Long64_t firstEv = 0, Long64_t maxEv = 1000000000,
                  TString simFile = "", int pdg = 0, double rCont = 145, double satAmp = 3500, double binLen = 15, double rExclAxis = 0, int pathMode = 1, double eomMaxBrho = 0.8, double braggMin = 2.0,
                  TString variants = "", int col = 1)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   std::map<int, double> gain;
   {
      std::ifstream fm(padMap.Length() ? padMap.Data() : "/dev/null");
      std::string line;
      while (std::getline(fm, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         int p; double g; long n; double x, y, gl;
         if (!(ss >> p >> g >> n)) continue;
         if (col == 4) { if (ss >> x >> y >> gl) gain[p] = gl; } else gain[p] = g;
      }
   }
   printf("pad map: %zu pads\n", gain.size());
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   std::vector<AtTools::AtHelix3D> hv;
   for (auto *o : *variants.Tokenize(";")) {
      hv.emplace_back();
      for (auto *kv : *((TObjString *)o)->GetString().Tokenize(",")) {
         TString k = ((TObjString *)kv)->GetString(), key = k(0, k.Index("=")), val = k(k.Index("=") + 1, k.Length());
         if (key == "hits") hv.back().SetHitRefit(val.Atoi() > 0, val.Atoi() > 1);
         else if (key == "cap") hv.back().SetEndWindowCap(val.Atof());
         else if (key == "axis") hv.back().SetAxisConstraint(val.Atof());
         else if (key == "guard") hv.back().SetRadiusGuard(val.Atof());
      }
   }
   std::vector<double> sv(hv.size()), bv(hv.size());
   std::vector<int> vv(hv.size());
   eom::Table tab[3];
   eom::makeTables(tab);
   namespace HX = AtTools::Helix3D;
   HX::EndOpt eo;
   eo.hitRefit = true; eo.hitRefitZ = true;
   TFile *fs = simFile.Length() ? TFile::Open(simFile) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   TClonesArray *mc = nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   TFile fo(out, "RECREATE");
   TTree o("w", "PID raw vs pad-gain corrected");
   int nh, va, vp, v3, ok1;
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s0, b0, s1, b1;
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s0", &s0); o.Branch("b0", &b0); o.Branch("s1", &s1);
   o.Branch("b1", &b1); o.Branch("ok1", &ok1);
   double hs, hb, hth, tbr, lseg, ta, bA, ef, em;
   int en, eit;
   o.Branch("hs", &hs); o.Branch("hb", &hb); o.Branch("hth", &hth); o.Branch("tbr", &tbr); o.Branch("en", &en);
   o.Branch("lseg", &lseg); o.Branch("ta", &ta); o.Branch("bA", &bA); o.Branch("ef", &ef); o.Branch("em", &em);
   o.Branch("eit", &eit);
   double bm, bke, bchi, bsat, bmc, brr;
   int bcont, bndf;
   o.Branch("bcont", &bcont); o.Branch("bm", &bm); o.Branch("bke", &bke); o.Branch("bchi", &bchi); o.Branch("bndf", &bndf);
   o.Branch("bsat", &bsat); o.Branch("bmc", &bmc); o.Branch("brr", &brr);
   double ef2, em2, efc, emc;
   int eit2, ksel;
   o.Branch("ef2", &ef2); o.Branch("em2", &em2); o.Branch("eit2", &eit2); o.Branch("efc", &efc); o.Branch("emc", &emc);
   o.Branch("ksel", &ksel);
   std::vector<double> vps, vpq;
   std::vector<int> vpu;
   double tke = -1;
   o.Branch("ps", &vps); o.Branch("pq", &vpq); o.Branch("pu", &vpu); o.Branch("tke", &tke);
   for (size_t k = 0; k < hv.size(); ++k) {
      o.Branch(Form("sv%zu", k), &sv[k]); o.Branch(Form("bv%zu", k), &bv[k]); o.Branch(Form("vv%zu", k), &vv[k]);
   }
   for (Long64_t ev = firstEv; ev < std::min(firstEv + maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      tbr = -1;
      if (ts) {
         int nl = 0;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (!m || m->GetPdgCode() != pdg) continue;
            tbr = std::sqrt(m->GetPx() * m->GetPx() + m->GetPy() * m->GetPy() + m->GetPz() * m->GetPz()) / 0.299792458;
            { double pm = tbr * 299.792458, mm = pdg == 2212 ? 938.272 : pdg == 1000010020 ? 1875.613 : 2808.921; tke = std::sqrt(pm * pm + mm * mm) - mm; }
            break;
         }
         if (tbr <= 0) continue;
      }
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         nh = H.size();
         if (nh < 30) continue;
         // pad-corrected copy of the track
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
         s0 = r0.valid ? r0.sqrtdEdx : -1; b0 = r0.valid ? r0.brho : -1;
         s1 = r1.valid ? r1.sqrtdEdx : -1; b1 = r1.valid ? r1.brho : -1;
         for (size_t k = 0; k < hv.size(); ++k) {
            auto rv = hv[k].Estimate(tc, sr.valid ? sr.direction : -1);
            vv[k] = rv.valid; sv[k] = rv.valid ? rv.sqrtdEdx : -1; bv[k] = rv.valid ? rv.brho : -1;
         }
         // hit-refit vertex-end helix (as AtHelix3D, vertex end from the 2-of-3 vote), arc angle, species-blind EOM
         hs = hb = hth = lseg = ta = bA = ef = em = -1; en = 0; eit = 0;
         bcont = 0; bm = bke = bchi = bsat = bmc = brr = -1; bndf = 0;
         ef2 = em2 = efc = emc = -1; eit2 = 0; ksel = 0;
         vps.clear(); vpq.clear(); vpu.clear();
         {
            std::vector<HX::P> h;
            for (auto &x : tc.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
            auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
            if (r.ok) {
               const double zA = h[r.order.front()].z, zB = h[r.order.back()].z;
               int spyEnd = -1;
               if (sr.valid && (sr.direction == 0 || sr.direction == 1)) spyEnd = (sr.direction == 0) == (zA < zB) ? 0 : 1;
               auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
               auto &E = r.e[V.end];
               double bq, sq, th, path;
               if (HX::firstArcPID(r, V.end, 152.0, bq, sq, th, path)) { hs = sq; hb = bq; hth = th; lseg = path; }
               const int mC = r.cen.size();
               std::vector<double> hx, hy, hz;
               if (E.ok && E.nUsed >= 1) {
                  const double lim = (V.end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[mC - E.nUsed]) + 2.0;
                  std::vector<std::pair<double, int>> w;
                  for (size_t i = 0; i < r.order.size(); ++i) {
                     double sf = V.end == 0 ? r.s[i] : r.length - r.s[i];
                     if (sf <= lim) w.push_back({sf, r.order[i]});
                  }
                  std::sort(w.begin(), w.end());
                  for (auto &pr : w) { hx.push_back(h[pr.second].x); hy.push_back(h[pr.second].y); hz.push_back(1000.0 - h[pr.second].z); }
               }
               en = hx.size();
               // arc-length angle: z vs R*|phi| round the end circle, first 50 % of the window hits
               if (en >= 12) {
                  int n = std::max(6, en / 2);
                  double prev = std::atan2(hy[0] - E.yc, hx[0] - E.xc), ph = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
                  for (int i = 0; i < n; ++i) {
                     double a = std::atan2(hy[i] - E.yc, hx[i] - E.xc), dph = a - prev;
                     while (dph > M_PI) dph -= 2 * M_PI;
                     while (dph < -M_PI) dph += 2 * M_PI;
                     if (i) ph += dph;
                     prev = a;
                     double sl = E.R * std::abs(ph);
                     sx += sl; sy += hz[i]; sxx += sl * sl; sxy += sl * hz[i];
                  }
                  double den = n * sxx - sx * sx;
                  if (std::abs(den) > 1e-9) { ta = std::atan2(1.0, (n * sxy - sx * sy) / den); bA = 0.00285 * E.R / std::sin(ta); }
               }
               // BRAGG FIT for tracks that stop on the small pads
               const auto &cl = r.cen[V.end == 0 ? mC - 1 : 0]; // far end
               if (kcal > 0 && mC >= 6 && std::hypot(cl.x, cl.y) < rCont && cl.z > 5 && cl.z < 995) {
                  bcont = 1;
                  // centroid order from the vertex end; saturation flags by re-binning exactly as the core does
                  std::vector<char> sat(mC, 0);
                  {
                     const auto &raw = tr.GetHitArray();
                     int c = 0;
                     for (size_t i = 0; i < r.order.size() && c < mC;) {
                        size_t j = i;
                        while (j < r.order.size() && r.s[j] - r.s[i] < 4.0) {
                           if (raw[r.order[j]]->GetCharge() > satAmp) sat[c] = 1;
                           ++j;
                        }
                        ++c; i = j;
                     }
                  }
                  auto at = [&](int jj) { return V.end == 0 ? jj : mC - 1 - jj; };
                  // smoothed (5-point) centroid positions; path = xy length / sin(theta_arc) (pathMode 1) so the
                  // centroid z jitter cannot add length -- the 3D polyline made slow ~90 deg tracks up to 2x longer
                  // than their range (elastic deuterons, 2026-11-01); pathMode 0 = old 3D polyline (3-point)
                  std::vector<double> px(mC), py(mC), pz(mC);
                  const int hw = pathMode == 1 ? 2 : 1;
                  const double sinT = (pathMode == 1 && ta > 0) ? std::sin(ta) : -1;
                  for (int jj = 0; jj < mC; ++jj) {
                     double sx = 0, sy = 0, sz = 0; int nn = 0;
                     for (int d = -hw; d <= hw; ++d) {
                        int kk = jj + d;
                        if (kk < 0 || kk >= mC) continue;
                        sx += r.cen[at(kk)].x; sy += r.cen[at(kk)].y; sz += r.cen[at(kk)].z; ++nn;
                     }
                     px[jj] = sx / nn; py[jj] = sy / nn; pz[jj] = sz / nn;
                  }
                  bragg::Profile P;
                  P.k = kcal; P.tp = &tab[0];
                  double L = 0, qsum = 0;
                  for (int jj = 0; jj < mC; ++jj) {
                     if (jj) {
                        double dxy = std::hypot(px[jj] - px[jj - 1], py[jj] - py[jj - 1]);
                        L += sinT > 0.05 ? dxy / sinT : std::sqrt(dxy * dxy + (pz[jj] - pz[jj - 1]) * (pz[jj] - pz[jj - 1]));
                     }
                     const auto &cc = r.cen[at(jj)];
                     P.s.push_back(L); P.q.push_back(cc.q); qsum += cc.q;
                     const double rc = std::hypot(cc.x, cc.y);
                     P.use.push_back(rc < rExclAxis ? 2 : (!sat[at(jj)] && rc <= 152.0 ? 1 : 0));
                  }
                  const int n = P.s.size();
                  for (int jj = 0; jj < n; ++jj) {
                     P.a.push_back(jj == 0 ? 0.0 : 0.5 * (P.s[jj - 1] + P.s[jj]));
                     P.b.push_back(jj == n - 1 ? P.s[jj] + 0.5 * (P.s[jj] - P.s[jj - 1]) : 0.5 * (P.s[jj] + P.s[jj + 1]));
                  }
                  P = bragg::rebin(P, binLen, brr);
                  vps = P.s; vpq = P.q;
                  for (char u : P.use) vpu.push_back(u);
                  int nsat = 0;
                  for (int jj = 0; jj < mC; ++jj) nsat += sat[jj];
                  bsat = double(nsat) / mC;
                  double ke0 = kcal * qsum, m0 = 1875.6;
                  if (bA > 0) {
                     double p = bA * 299.792458;
                     double mc = (p * p - ke0 * ke0) / (2 * ke0);
                     bmc = mc;
                     if (mc > 500 && mc < 6000) m0 = mc;
                  }
                  double mf = m0, kf = ke0, c2;
                  int nd;
                  if (bragg::fit(P, mf, kf, c2, nd)) { bm = mf; bke = kf; bchi = c2; bndf = nd; }
               }
               // species-blind EOM, two paths
               const bool okEom = kcal > 0 && bA > 0 && en >= 10 && ta > 0 && !E.axisConstrained && !E.guarded;
               if (okEom && hs > 0 && lseg > 5 && bA < eomMaxBrho) { // through-going path: mass from the mean dE/dx
                  const double Smeas = kcal * hs * hs;
                  double p0 = bA * 299.792458, m = eom::solveMass(tab[0], p0, lseg, Smeas);
                  for (eit = 1; eit <= 3; ++eit) {
                     double fm, tho;
                     double pn = eom::fitM(hx, hy, hz, E.xc, E.yc, E.R, E.sense, ta, m, tab[0], 0, fm, tho);
                     if (pn <= 0) { p0 = -1; break; }
                     p0 = pn;
                     double mn = eom::solveMass(tab[0], p0, lseg, Smeas);
                     bool conv = std::abs(mn / m - 1) < 0.005;
                     m = mn;
                     if (conv) break;
                  }
                  if (p0 > 0) { ef = p0 / 299.792458; em = m; }
               }
               if (okEom && bm > 0 && bke > 0.02) { // stopping path (every contained track with a Bragg fit): mass from the Bragg energy
                  const double KE = bke;
                  auto massOf = [&](double p) { return (p * p - KE * KE) / (2 * KE); };
                  double p0 = bA * 299.792458, m = massOf(p0);
                  if (m > 300 && m < 8000) {
                     for (eit2 = 1; eit2 <= 3; ++eit2) {
                        double fm, tho;
                        double pn = eom::fitM(hx, hy, hz, E.xc, E.yc, E.R, E.sense, ta, m, tab[0], 0, fm, tho);
                        if (pn <= 0) { p0 = -1; break; }
                        p0 = pn;
                        double mn = massOf(p0);
                        if (mn < 300 || mn > 8000) { p0 = -1; break; }
                        bool conv = std::abs(mn / m - 1) < 0.005;
                        m = mn;
                        if (conv) break;
                     }
                     if (p0 > 0) { ef2 = p0 / 299.792458; em2 = m; }
                  }
               }
               // combined: the stopping path if the track shows a Bragg rise OR the dE/dx solver failed (no mass, or railed
               // at its upper limit -- a track that stops inside the PID segment, 2026-11-01 truth test); else through-going
               const bool emOk = em > 0 && em < 7900;
               if (em2 > 0 && (brr > braggMin || !emOk)) { emc = em2; efc = ef2; ksel = 2; }
               else if (emOk) { emc = em; efc = ef; ksel = 1; }
            }
         }
         b3 = b0; sq3h = s0; th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
