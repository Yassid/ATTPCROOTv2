/// @file pid_dedxeom_dp.C
/// @brief SPECIES-BLIND equation-of-motion Brho (2026-11-01): the energy loss along the trajectory follows the
///        track's OWN measured dE/dx, not a mass hypothesis. For Z = 1 the stopping power depends on velocity
///        only, so the measured mean dE/dx over the PID segment (k * charge/mm, k [MeV per charge unit] calibrated
///        per run on the deuteron band by dedxeom_calib_dp.C) plus the momentum fixes the mass; iterate:
///          m <- solveMass(p0, Lseg, k*hs^2);  p0 <- planar EOM fit (velocity-scaled proton CATIMA table, mass m,
///          theta fixed to the arc-length angle).  3 iterations (stop at |dm/m| < 0.5 %).
///        Branches as pid_padgain_test_dp.C plus: hs (pad-corrected sqrt charge/mm of the AtHelix3D hit-refit PID
///        segment), lseg (its helix path, mm), hb/hth (hit-refit Brho, theta), ta (arc-length angle: z vs R*phi,
///        first 50 % of the window), bA (hit-refit R with ta), and with k > 0: ef (EOM Brho, T m), em (fitted mass,
///        MeV), eit (iterations). k <= 0: calibration pass (no fits).
///   root -l -b -q 'pid_dedxeom_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0026_multifit_reco.root","pid/padgain_d2_v1.txt","out.root",0)'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
#include "eom_brho_dp.h"
void pid_dedxeom_dp(TString in, TString padMap, TString out, double kcal, Long64_t firstEv = 0, Long64_t maxEv = 1000000000,
                    TString simFile = "", int pdg = 0, double eomMaxBrho = 0.8, TString variants = "", int col = 1)
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
               // species-blind EOM: mass from the measured dE/dx, iterated with the trajectory fit
               if (kcal > 0 && hs > 0 && lseg > 5 && bA > 0 && bA < eomMaxBrho && en >= 10 && !E.axisConstrained && !E.guarded) {
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
            }
         }
         b3 = b0; sq3h = s0; th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
