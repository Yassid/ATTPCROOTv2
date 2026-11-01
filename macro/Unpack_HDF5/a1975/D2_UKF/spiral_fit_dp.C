/// @file spiral_fit_dp.C
/// @brief SHRINKING-RADIUS (spiral) vertex-end fit, prototype on simulation (2026-11-01). A decelerating track's xy
///        projection tightens: instead of one circle, fit r(phi) = R0 + k1*phi (+ k2*phi^2) round a free centre, phi =
///        unwrapped |phase| from the vertex-end hit, so R0 is the radius AT THE VERTEX. Species-blind (k free).
///        Polar angle: z vs arc length s = integral r(phi) dphi, first 50 % of the window hits. Brho = 2.85e-3 R0 / sin(theta).
///        Windows: w0 = the AtHelix3D hit-refit window (as bA), w1 = 2x that path, w2 = from the vertex end up to one turn.
///        Writes per track: truth tbr, Spyral br/pol, bA/ta (reproduced), and for each window/model the spiral Brho.
///   root -l -b -q 'spiral_fit_dp.C("/mnt/f/a1975_C16_dt_sim/gs_s3001_reco.root","/mnt/f/a1975_C16_dt_sim/gs_s3001_sim.root","out.root",1000010030)'
///   data: spiral_fit_dp.C("<reco>","","out.root",0,maxEv,"pid/padgain_d2_v1.txt",firstEv)
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
namespace spf {
std::vector<double> X, Y;
int nPar = 4;
void fcn(int &, double *, double &f, double *p, int)
{
   // p: xc, yc, R0, k1, k2
   double prev = std::atan2(Y[0] - p[1], X[0] - p[0]), ph = 0, c2 = 0;
   for (size_t i = 0; i < X.size(); ++i) {
      double a = std::atan2(Y[i] - p[1], X[i] - p[0]), d = a - prev;
      while (d > M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      if (i) ph += d;
      prev = a;
      double aph = std::abs(ph);
      double rm = p[2] + p[3] * aph + (nPar > 4 ? p[4] * aph * aph : 0);
      double res = std::hypot(X[i] - p[0], Y[i] - p[1]) - rm;
      c2 += res * res;
   }
   f = c2;
}
/// fit; returns R0 (and k1, k2, centre); false on failure
bool fit(const std::vector<double> &x, const std::vector<double> &y, double xc0, double yc0, double R00, int np, double *out)
{
   X = x; Y = y; nPar = np;
   TMinuit mn(np);
   mn.SetPrintLevel(-1);
   mn.SetFCN(fcn);
   int ierr = 0;
   double arg[2] = {1, 0};
   mn.mnexcm("SET ERR", arg, 1, ierr);
   mn.mnparm(0, "xc", xc0, 1, 0, 0, ierr);
   mn.mnparm(1, "yc", yc0, 1, 0, 0, ierr);
   mn.mnparm(2, "R0", R00, 1, 1, 5000, ierr);
   mn.mnparm(3, "k1", 0, 0.5, 0, 0, ierr);
   if (np > 4) mn.mnparm(4, "k2", 0, 0.1, 0, 0, ierr);
   arg[0] = 5000; arg[1] = 0.01;
   mn.mnexcm("MIGRAD", arg, 2, ierr);
   double e;
   for (int i = 0; i < np; ++i) mn.GetParameter(i, out[i], e);
   if (np == 4) out[4] = 0;
   double fmin, edm, errdef;
   int nv, nx, ist;
   mn.mnstat(fmin, edm, errdef, nv, nx, ist);
   return ist >= 1 && out[2] > 1;
}
/// polar angle from z vs arc length along the fitted spiral, first frac of the hits
double theta(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, const double *p, double frac = 0.5)
{
   int n = std::min<int>(x.size(), std::max(6, (int)(frac * x.size())));
   double prev = std::atan2(y[0] - p[1], x[0] - p[0]), ph = 0, s = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
   for (int i = 0; i < n; ++i) {
      double a = std::atan2(y[i] - p[1], x[i] - p[0]), d = a - prev;
      while (d > M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      if (i) ph += d;
      prev = a;
      // arc length from the NET phase (summing |dphi| would inflate s with the hit-order jitter)
      double aph = std::abs(ph);
      s = p[2] * aph + 0.5 * p[3] * aph * aph + p[4] * aph * aph * aph / 3.0;
      sx += s; sy += z[i]; sxx += s * s; sxy += s * z[i];
   }
   double den = n * sxx - sx * sx;
   if (std::abs(den) < 1e-9) return -1;
   return std::atan2(1.0, (n * sxy - sx * sy) / den);
}
} // namespace spf
void spiral_fit_dp(TString in, TString simFile, TString out, int pdg, Long64_t maxEv = 1000000000, TString padMap = "",
                   Long64_t firstEv = 0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   namespace HX = AtTools::Helix3D;
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   // data: no simFile -> every track with >= 30 hits of every event, in the order pid_sbpid_dp.C writes them
   // (friend-able with its pv_<run>.root); pad-gain map applied to the charges as there
   TFile *fs = simFile.Length() ? TFile::Open(simFile) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   TClonesArray *mc = nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   std::map<int, double> gain;
   {
      std::ifstream fm(padMap.Length() ? padMap.Data() : "/dev/null");
      std::string line;
      while (std::getline(fm, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         int pp; double g; long n;
         if (ss >> pp >> g >> n) gain[pp] = g;
      }
   }
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   HX::EndOpt eo;
   eo.hitRefit = true; eo.hitRefitZ = true;
   TFile fo(out, "RECREATE");
   TTree o("w", "spiral vertex-end fit");
   double tbr, br, pol, bA, ta, R, lwin;
   int nw[3];
   double sb[3][2], sR0[3][2], sk1[3][2], sth[3][2];
   o.Branch("tbr", &tbr); o.Branch("br", &br); o.Branch("pol", &pol); o.Branch("bA", &bA); o.Branch("ta", &ta);
   o.Branch("R", &R); o.Branch("lwin", &lwin);
   for (int w = 0; w < 3; ++w) {
      o.Branch(Form("n%d", w), &nw[w]);
      for (int m = 0; m < 2; ++m) {
         o.Branch(Form("b%d%d", w, m), &sb[w][m]); o.Branch(Form("R%d%d", w, m), &sR0[w][m]);
         o.Branch(Form("k%d%d", w, m), &sk1[w][m]); o.Branch(Form("t%d%d", w, m), &sth[w][m]);
      }
   }
   long nDone = 0;
   for (Long64_t ev = firstEv; ev < std::min(firstEv + maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto &cands = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      tbr = -1;
      if (ts) {
         int nl = 0;
         for (auto &tr : cands) nl += tr.GetHitArray().size() >= 30;
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
      for (auto &tr0 : cands) {
         if (tr0.GetHitArray().size() < 30) continue;
         auto sr = spy.Estimate(tr0);
         AtTrack tr(tr0);
         for (auto &hh : tr.GetHitArray()) {
            auto it = gain.find(hh->GetPadNum());
            if (it != gain.end() && it->second > 0) hh->SetCharge(hh->GetCharge() / it->second);
         }
         br = sr.valid ? sr.brho : -1; pol = sr.valid ? sr.polar * TMath::RadToDeg() : -999;
         bA = ta = R = lwin = -1;
         for (int w = 0; w < 3; ++w) { nw[w] = 0; for (int m = 0; m < 2; ++m) sb[w][m] = sR0[w][m] = sk1[w][m] = sth[w][m] = -1; }
         std::vector<HX::P> h;
         for (auto &x : tr.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
         if (!r.ok || r.length / std::max(1, r.nBackbone) >= 6) { o.Fill(); continue; }
         const double zA = h[r.order.front()].z, zB = h[r.order.back()].z;
         int spyEnd = -1;
         if (sr.valid && (sr.direction == 0 || sr.direction == 1)) spyEnd = (sr.direction == 0) == (zA < zB) ? 0 : 1;
         auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
         auto &E = r.e[V.end];
         if (!E.ok || E.nUsed < 1) { o.Fill(); continue; }
         const int mC = r.cen.size();
         const double lim = (V.end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[mC - E.nUsed]) + 2.0;
         lwin = lim; R = E.R;
         // all hits ordered from the vertex end
         std::vector<std::pair<double, int>> all;
         for (size_t i = 0; i < r.order.size(); ++i) all.push_back({V.end == 0 ? r.s[i] : r.length - r.s[i], r.order[i]});
         std::sort(all.begin(), all.end());
         for (int w = 0; w < 3; ++w) {
            std::vector<double> x, y, z;
            double prev = 0, ph = 0;
            for (size_t i = 0; i < all.size(); ++i) {
               const auto &hp = h[all[i].second];
               if (w == 0 && all[i].first > lim) break;
               if (w == 1 && all[i].first > 2 * lim) break;
               if (w == 2) { // up to one turn round the end circle
                  double a = std::atan2(hp.y - E.yc, hp.x - E.xc);
                  if (i) { double d = a - prev; while (d > M_PI) d -= 2 * M_PI; while (d < -M_PI) d += 2 * M_PI; ph += d; }
                  prev = a;
                  if (std::abs(ph) > 2 * M_PI) break;
               }
               x.push_back(hp.x); y.push_back(hp.y); z.push_back(1000.0 - hp.z);
            }
            nw[w] = x.size();
            if (x.size() < 12) continue;
            if (w == 0) {
               // reproduce bA: circle + arc angle
               double pc[5] = {E.xc, E.yc, E.R, 0, 0};
               ta = spf::theta(x, y, z, pc);
               if (ta > 0) bA = 0.00285 * E.R / std::sin(ta);
            }
            for (int m = 0; m < 2; ++m) {
               double p[5];
               if (!spf::fit(x, y, E.xc, E.yc, E.R, m == 0 ? 4 : 5, p)) continue;
               double th = spf::theta(x, y, z, p);
               if (th <= 0) continue;
               sR0[w][m] = p[2]; sk1[w][m] = p[3]; sth[w][m] = th; sb[w][m] = 0.00285 * p[2] / std::sin(th);
            }
         }
         o.Fill();
         if (++nDone % 1000 == 0) printf("%ld tracks\n", nDone);
      }
   }
   fo.cd();
   o.Write();
   printf("done: %ld tracks\n", nDone);
}
