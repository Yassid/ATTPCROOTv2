/// @file pid_eom_dp.C
/// @brief Equation-of-motion Brho (eom_brho_dp.h) on data or simulation, same tracks and branches as
///        pid_variant_test_dp.C (sqa/br, sqp, s0/b0, s1/b1, reference v3 b3 sq3h th3 sp3, variants sv/bv/vv) plus,
///        from the AtHelix3D hit-refit (circle + z line) vertex-end helix: hs (its sqrt dE/dx, pad-corrected),
///        hb (its Brho), hth (theta) and the EOM fits for hypotheses k = 0 p, 1 d, 2 t:
///          e<m>b<k> (Brho [T m]), e<m>c<k> (chi2), e<m>t<k> (theta, mode 1) with m = 0 planar, 1 3D; en (hits fitted).
///        EOM only for tracks with hit-refit Brho < eomMaxBrho (no deceleration to model above).
///        Simulation: simFile + pdg add tbr (truth Brho); only events with exactly one >= 30-hit track are kept.
///   root -l -b -q 'pid_eom_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2_v1.txt","out.root","hits=2",0,5000)'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
#include "eom_brho_dp.h"
void pid_eom_dp(TString in, TString padMap, TString out, TString variants, Long64_t firstEv = 0, Long64_t maxEv = 1000000000,
                TString simFile = "", int pdg = 0, double eomMaxBrho = 1.2, int col = 1)
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
   double hs, hb, hth, tbr, eb[2][3], ec[2][3], et[2][3];
   int en;
   o.Branch("hs", &hs); o.Branch("hb", &hb); o.Branch("hth", &hth); o.Branch("tbr", &tbr); o.Branch("en", &en);
   for (int m = 0; m < 2; ++m)
      for (int k = 0; k < 3; ++k) {
         o.Branch(Form("e%db%d", m, k), &eb[m][k]); o.Branch(Form("e%dc%d", m, k), &ec[m][k]); o.Branch(Form("e%dt%d", m, k), &et[m][k]);
      }
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
         // hit-refit vertex-end helix (as AtHelix3D, vertex end from the 2-of-3 vote) and the EOM fits
         hs = hb = hth = -1; en = 0;
         for (int m = 0; m < 2; ++m) for (int k = 0; k < 3; ++k) eb[m][k] = ec[m][k] = et[m][k] = -1;
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
               if (HX::firstArcPID(r, V.end, 152.0, bq, sq, th, path)) { hs = sq; hb = bq; hth = th; }
               if (hb > 0 && hb < eomMaxBrho && !E.axisConstrained && !E.guarded) {
                  const int mC = r.cen.size();
                  const double lim = (V.end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[mC - E.nUsed]) + 2.0;
                  std::vector<std::pair<double, int>> w;
                  for (size_t i = 0; i < r.order.size(); ++i) {
                     double sf = V.end == 0 ? r.s[i] : r.length - r.s[i];
                     if (sf <= lim) w.push_back({sf, r.order[i]});
                  }
                  std::sort(w.begin(), w.end());
                  std::vector<double> hx, hy, hz;
                  for (auto &pr : w) { hx.push_back(h[pr.second].x); hy.push_back(h[pr.second].y); hz.push_back(1000.0 - h[pr.second].z); }
                  en = hx.size();
                  if (en >= 10)
                     for (int m = 0; m < 2; ++m)
                        for (int k = 0; k < 3; ++k) {
                           double fm, tho;
                           double p0 = eom::fit(hx, hy, hz, E.xc, E.yc, E.R, E.sense, E.theta, k, tab, m, fm, tho);
                           if (p0 > 0) { eb[m][k] = p0 / 299.792458; ec[m][k] = fm; et[m][k] = tho; }
                        }
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
