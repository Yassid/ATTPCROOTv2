/// @file pid_variants_dp.C
/// @brief Spyral PID plus several configurations of the 3D helix-path PID (helix3d.h), computed for every
///        pattern track in ONE pass, so the configurations are scored on identical tracks.
///
/// 3D PID = vertex end by the 2-of-3 majority (sense, Spyral, Bragg), Brho from that end's helix,
/// sqrt(dE/dx) = inlier charge over the first arc on the small pads / HELIX path (the estimator of
/// 2026-10-30). Variant k changes only the estimate() options listed in `cfg` below.
/// Branches: Spyral va br sqa pol; reference 3D (variant 0) v3 b3 sq3h th3 sp3; variants ok<k> b<k> s<k>.
/// The reference names match pid_width_dp.C so pid_proton_fixed.C / pid_dd_generic.C work unchanged.
///
///   root -l -b -q 'pid_variants_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","out.root")'
#include "../UKF/pid/helix3d.h"

namespace pv {
struct Cfg {
   const char *name;
   double binW, maxEdge;
   h3d::EndOpt eo;
};
// first-arc, small-pad, inlier charge over the helix path, from end `end` (same as pid_width_dp.C sq3h)
bool helixPID(const h3d::Result &r, int end, double &brho, double &sq, double &th)
{
   auto &E = r.e[end];
   if (!E.ok || E.brho <= 0) return false;
   const int m = r.cen.size();
   auto inl = [&](const h3d::P &p) { return std::abs(std::hypot(p.x - E.xc, p.y - E.yc) - E.R) < E.resCut; };
   int j0 = end == 0 ? 0 : m - 1, jFar = 0;
   for (int jj = 0; jj < m; ++jj) {
      int j = end == 0 ? jj : m - 1 - jj;
      if (inl(r.cen[j])) { j0 = j; break; }
   }
   double dFar = -1;
   for (int jj = 0; jj < m; ++jj) {
      int j = end == 0 ? jj : m - 1 - jj;
      if (!inl(r.cen[j])) continue;
      double d = std::hypot(r.cen[j].x - r.cen[j0].x, r.cen[j].y - r.cen[j0].y);
      if (d > dFar) { dFar = d; jFar = jj; }
   }
   double q = 0, phPrev = 0, phAcc = 0, phFirst = 0, phLast = 0;
   bool started = false;
   for (int jj = 0; jj <= jFar; ++jj) {
      int j = end == 0 ? jj : m - 1 - jj;
      auto &p = r.cen[j];
      double a = std::atan2(p.y - E.yc, p.x - E.xc);
      if (jj == 0) phPrev = a;
      else {
         double d = a - phPrev;
         while (d > M_PI) d -= 2 * M_PI;
         while (d < -M_PI) d += 2 * M_PI;
         phAcc += d; phPrev = a;
      }
      bool inner = std::hypot(p.x, p.y) <= 152.0;
      if (!inl(p)) continue;
      if (!inner && started) break;
      if (inner) {
         if (!started) { phFirst = phAcc; started = true; }
         phLast = phAcc;
         q += p.q;
      }
   }
   double path = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi) * std::abs(phLast - phFirst);
   if (!started || path < 5) return false;
   brho = E.brho; sq = std::sqrt(q / path); th = E.theta * TMath::RadToDeg();
   return true;
}
} // namespace pv

void pid_variants_dp(TString in, TString out, Long64_t maxEv = 1000000000)
{
   using namespace pv;
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   std::vector<Cfg> cfg;
   {
      h3d::EndOpt d; // the ADOPTED defaults (onC6, axis dca>10) -- binW/maxEdge passed explicitly below
      cfg.push_back({"adopted default", 4.0, -4.0, d});
   }
   const int K = cfg.size();
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   TFile fo(out, "RECREATE");
   TTree o("w", "Spyral + 3D PID variants");
   int nh, va, v3, ok[16];
   double br, sqa, pol, b3, sq3h, th3, sp3, b[16], s[16];
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("br", &br); o.Branch("sqa", &sqa); o.Branch("pol", &pol);
   o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h); o.Branch("th3", &th3); o.Branch("sp3", &sp3);
   for (int k = 0; k < K; ++k) {
      o.Branch(Form("ok%d", k), &ok[k]); o.Branch(Form("b%d", k), &b[k]); o.Branch(Form("s%d", k), &s[k]);
   }
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         nh = H.size();
         if (nh < 30) continue;
         std::vector<h3d::P> h;
         for (auto &x : H) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto sr = spy.Estimate(tr);
         va = sr.valid; br = sr.valid ? sr.brho : -1; sqa = sr.valid ? sr.sqrtdEdx : -1;
         pol = sr.valid ? sr.polar * TMath::RadToDeg() : -999;
         v3 = 0; b3 = sq3h = th3 = -1; sp3 = 1e9;
         for (int k = 0; k < K; ++k) {
            ok[k] = 0; b[k] = s[k] = -1;
            auto r = h3d::estimate(h, 2.85, -1, 1000, 5, cfg[k].binW, cfg[k].maxEdge, 8, 0.15, 0.5, cfg[k].eo);
            if (!r.ok) continue;
            int spyEnd = -1; // Spyral forward = vertex at the LOWEST z_digi; h is digi frame
            if (sr.valid) spyEnd = (sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1;
            auto V = h3d::majorityEnd(r, spyEnd);
            double bb, ss, th;
            if (!helixPID(r, V.end, bb, ss, th)) continue;
            ok[k] = 1; b[k] = bb; s[k] = ss;
            if (k == 0) { v3 = 1; b3 = bb; sq3h = ss; th3 = th; sp3 = r.length / r.nBackbone; }
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks;", out.Data(), o.GetEntries());
   for (int k = 0; k < K; ++k) printf("  v%d=%s", k, cfg[k].name);
   printf("\n");
}
