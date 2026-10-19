/// @file census_arc_dp.C -- PRE-Spyral census of 16C(d,p) pattern tracks, data vs g.s. sim (2026-10-12).
/// Question: the data keep ~0.3 of the sim's reconstructed forward-cm bound protons at vz 200-400 (exiting arcs).
/// effchain_bound_dp.C showed Spyral-valid -> final keeps 0.68-0.90 flat, so the loss is BEFORE that: either the
/// pattern track is missing, or AtSpyralPID rejects it (failCode is transient, so it is re-run here).
/// Per pattern track: an independent helix estimate from the hits (Kasa circle in xy + linear z vs unwrapped phase),
/// and AtSpyralPID::Estimate with B=2.85 and default settings (as AtGenfitter builds it).
/// Frames: hit z is the digi frame; vz_fit = 1000 - z_hit, theta_fit = 180 - theta_hit.
/// Usage: root -l -b -q 'census_arc_dp.C("in_reco.root","out.root",src,run)'   src 0 data, 1 sim
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz) { return 184.25 - 0.011030 * (vz - 500); }

struct Helix {
   bool ok{false};
   double xc, yc, R, s, thHit, zv, d0;
};
// Kasa circle + z(phi) line, hits ordered by z
Helix helix1(std::vector<std::array<double, 4>> &h)
{
   Helix H;
   int n = h.size();
   if (n < 5) return H;
   TMatrixD A(n, 3);
   TVectorD b(n);
   for (int i = 0; i < n; ++i) {
      A(i, 0) = h[i][0]; A(i, 1) = h[i][1]; A(i, 2) = 1;
      b(i) = h[i][0] * h[i][0] + h[i][1] * h[i][1];
   }
   TMatrixD At(TMatrixD::kTransposed, A);
   TMatrixD N = At * A;
   if (std::abs(N.Determinant()) < 1e-12) return H;
   TVectorD c = N.Invert() * (At * b);
   H.xc = c(0) / 2; H.yc = c(1) / 2;
   double r2 = c(2) + H.xc * H.xc + H.yc * H.yc;
   if (r2 <= 0) return H;
   H.R = std::sqrt(r2);
   std::vector<double> ph(n);
   for (int i = 0; i < n; ++i) {
      ph[i] = std::atan2(h[i][1] - H.yc, h[i][0] - H.xc);
      if (i) { while (ph[i] - ph[i - 1] > M_PI) ph[i] -= 2 * M_PI; while (ph[i] - ph[i - 1] < -M_PI) ph[i] += 2 * M_PI; }
   }
   double sp = 0, sz = 0, spp = 0, spz = 0;
   for (int i = 0; i < n; ++i) { sp += ph[i]; sz += h[i][2]; spp += ph[i] * ph[i]; spz += ph[i] * h[i][2]; }
   double den = n * spp - sp * sp;
   if (den <= 0) return H;
   H.s = (n * spz - sp * sz) / den;
   double z0 = (sz - H.s * sp) / n;
   H.thHit = std::atan2(H.R, std::abs(H.s)); // 0..90 deg w.r.t. +z_hit
   // vertex: circle point nearest the axis, phase taken on the branch nearest the zmin-end hit
   double pv = std::atan2(-H.yc, -H.xc);
   while (pv - ph[0] > M_PI) pv -= 2 * M_PI;
   while (pv - ph[0] < -M_PI) pv += 2 * M_PI;
   H.zv = z0 + H.s * pv;
   H.d0 = std::abs(std::hypot(H.xc, H.yc) - H.R);
   H.ok = true;
   return H;
}
// Two passes: the whole track, then only the first HALF TURN from the vertex end (the spiral shrinks as the
// particle slows, so a whole-track circle reads the momentum ~10 % low). vEnd: 0 = zmin end, 1 = zmax end.
Helix helix(std::vector<std::array<double, 4>> &h, int vEnd)
{
   std::sort(h.begin(), h.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
   std::vector<std::array<double, 4>> w = h;
   if (vEnd) std::reverse(w.begin(), w.end());
   Helix H = helix1(w);
   if (!H.ok) return H;
   double p0 = std::atan2(w[0][1] - H.yc, w[0][0] - H.xc), prev = p0, acc = 0;
   std::vector<std::array<double, 4>> f;
   for (auto &x : w) {
      double p = std::atan2(x[1] - H.yc, x[0] - H.xc), d = p - prev;
      while (d > M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      acc += d; prev = p;
      if (std::abs(acc) > M_PI) break;
      f.push_back(x);
   }
   if (f.size() >= 8 && f.size() < w.size()) {
      Helix H2 = helix1(f);
      if (H2.ok) H = H2;
   }
   return H;
}

void census_arc_dp(TString in, TString out, int src, int run)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   TFile *f = TFile::Open(in);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pa);
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   TFile fo(out, "RECREATE");
   TTree o("c", "pattern-track census");
   int S = src, RUN = run, ev, tid, nh, ncl, sv, sf, dirB;
   double zmin, zmax, rzmin, rzmax, qs, R, d0, thl, vz, ke, ex, tcm, sth, svz, ske, sex, stcm, ssq, sbr;
   o.Branch("src", &S); o.Branch("run", &RUN); o.Branch("ev", &ev); o.Branch("tid", &tid); o.Branch("nh", &nh);
   o.Branch("ncl", &ncl); o.Branch("zmin", &zmin); o.Branch("zmax", &zmax); o.Branch("rzmin", &rzmin);
   o.Branch("rzmax", &rzmax); o.Branch("qs", &qs); o.Branch("R", &R); o.Branch("d0", &d0); o.Branch("dirB", &dirB);
   o.Branch("thl", &thl); o.Branch("vz", &vz); o.Branch("ke", &ke); o.Branch("ex", &ex); o.Branch("tcm", &tcm);
   o.Branch("sv", &sv); o.Branch("sf", &sf); o.Branch("sth", &sth); o.Branch("svz", &svz); o.Branch("ske", &ske);
   o.Branch("sex", &sex); o.Branch("stcm", &stcm); o.Branch("ssq", &ssq); o.Branch("sbr", &sbr);
   const double mp = 938.272;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      if (t->GetEntry(i) <= 0 || !pa || pa->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)pa->At(0);
      if (!pe) continue;
      ev = i;
      for (auto &tr : pe->GetTrackCand()) {
         auto &hits = tr.GetHitArray();
         nh = hits.size();
         if (nh < 5) continue;
         tid = tr.GetTrackID();
         ncl = tr.GetHitClusterArray() ? tr.GetHitClusterArray()->size() : -1;
         std::vector<std::array<double, 4>> h;
         qs = 0;
         for (auto &ht : hits) { auto p = ht->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), ht->GetCharge()}); qs += ht->GetCharge(); }
         std::sort(h.begin(), h.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
         zmin = h.front()[2]; zmax = h.back()[2];
         rzmin = std::hypot(h.front()[0], h.front()[1]); rzmax = std::hypot(h.back()[0], h.back()[1]);
         dirB = rzmin < rzmax; // vertex end at low z_hit = moving to +z_hit = backward lab
         auto H = helix(h, dirB ? 0 : 1);
         R = H.ok ? H.R : -1; d0 = H.ok ? H.d0 : -1;
         thl = vz = ke = ex = tcm = -999;
         if (H.ok) {
            thl = dirB ? 180 - H.thHit * TMath::RadToDeg() : H.thHit * TMath::RadToDeg();
            vz = 1000 - H.zv;
            double p = 0.299792458 * 2.85 * H.R / std::sin(H.thHit);
            ke = std::sqrt(p * p + mp * mp) - mp;
            if (!kine2b(Eb(vz), thl * TMath::DegToRad(), ke, ex, tcm)) ex = tcm = -999;
         }
         auto r = spy.Estimate(tr);
         sv = r.valid; sf = r.failCode;
         sth = svz = ske = sex = stcm = ssq = sbr = -999;
         if (r.valid) {
            sth = 180 - r.polar * TMath::RadToDeg(); svz = 1000 - r.vertex.Z();
            double p = 299.792458 * r.brho;
            ske = (std::sqrt(p * p + mp * mp) - mp) / 1.023; ssq = r.sqrtdEdx; sbr = r.brho;
            if (!kine2b(Eb(svz), sth * TMath::DegToRad(), ske, sex, stcm)) sex = stcm = -999;
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("done %s -> %s\n", in.Data(), out.Data());
}
