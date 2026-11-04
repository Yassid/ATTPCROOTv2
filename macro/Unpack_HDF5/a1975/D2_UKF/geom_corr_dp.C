#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
/// @file geom_corr_dp.C
/// Which single variable orders the (d,p) loss band? For each variable X (per true proton, from geom_vars_dp.C or the sim census), fits the (vz, theta_cm) cells with
/// data = F(theta_cm) x sum_k g(k) sim(cell, k), k = quantile bin of X, Poisson likelihood. F free per angle, so any dsigma shape is allowed; g is the loss. Deviance/ndf ~ 1 means X explains the band.
///   root -l -b -q geom_corr_dp.C
namespace gc { using namespace exsel;
const int NT = 5, NZ = 7, NC = NT * NZ, KQ = 5; const double TB[NT + 1] = {10, 20, 25, 30, 35, 40}, ZB[NZ + 1] = {200, 300, 400, 500, 600, 700, 800, 900};
int tb(double t) { for (int i = 0; i < NT; ++i) if (t >= TB[i] && t < TB[i + 1]) return i; return -1; }
int zb(double z) { for (int i = 0; i < NZ; ++i) if (z >= ZB[i] && z < ZB[i + 1]) return i; return -1; }
double EbB(double vz) { return 183.873 - 0.011099 * vz; }
double dev(const double *d, const double *m) { double D = 0; for (int c = 0; c < NC; ++c) { D += 2 * (m[c] - d[c]); if (d[c] > 0) D += 2 * d[c] * log(d[c] / m[c]); else if (m[c] < 1e-20) continue; } return D; }
struct Res { double D, F[NT], g[KQ], q[KQ + 1], mu[NC]; };
/// multiplicative fit; s[c][k] sim tracks in cell c, quantile k
Res fit(const double *d, double s[NC][KQ], int nk) { Res R; double F[NT], g[KQ]; for (int j = 0; j < NT; ++j) F[j] = 0.2; for (int k = 0; k < KQ; ++k) g[k] = 1;
   for (int it = 0; it < 4000; ++it) { double mu[NC]; for (int c = 0; c < NC; ++c) { double a = 1e-9; for (int k = 0; k < nk; ++k) a += g[k] * s[c][k]; mu[c] = std::max(F[c % NT] * a, 1e-30); }
      double num[KQ] = {0}, den[KQ] = {0}; for (int c = 0; c < NC; ++c) for (int k = 0; k < nk; ++k) { den[k] += F[c % NT] * s[c][k]; num[k] += d[c] * F[c % NT] * g[k] * s[c][k] / mu[c]; }
      for (int k = 0; k < nk; ++k) if (den[k] > 0) g[k] = std::max(num[k] / den[k], 1e-6);
      double top = 0; for (int k = 0; k < nk; ++k) top = std::max(top, g[k]); for (int k = 0; k < nk; ++k) g[k] /= top;
      for (int j = 0; j < NT; ++j) { double n = 0, dn = 0; for (int c = 0; c < NC; ++c) if (c % NT == j) { double a = 1e-9; for (int k = 0; k < nk; ++k) a += g[k] * s[c][k]; n += d[c]; dn += a; } F[j] = n / dn; } }
   for (int c = 0; c < NC; ++c) { double a = 1e-9; for (int k = 0; k < nk; ++k) a += g[k] * s[c][k]; R.mu[c] = std::max(F[c % NT] * a, 1e-30); } R.D = dev(d, R.mu); for (int j = 0; j < NT; ++j) R.F[j] = F[j]; for (int k = 0; k < nk; ++k) R.g[k] = g[k]; return R; }
struct T { int cell; std::vector<double> x; };
const int NV = 44;
bool masked(int c, int mask) { int a = c % NT, z = c / NT; return (mask == 1 && a == 0) || (mask == 2 && a == 0 && z == NZ - 1) || (mask == 3 && (a == 0 || a == 4)); }
/// per-cell picture: data/sim, and the sim medians of the leading variables
void cellmap(const std::vector<T> &tr, const double *d, const std::vector<std::string> &N) {
   double s[NC] = {0}, dn = 0, sn = 0; for (auto &t : tr) s[t.cell]++; for (int c = 0; c < NC; ++c) if (c % NT >= 3) { dn += d[c]; sn += s[c]; }
   const char *V[6] = {"path", "zext", "ring100-150", "ring60-100", "turns", "rmax"}; std::vector<int> iv; for (auto n : V) iv.push_back(std::find(N.begin(), N.end(), n) - N.begin());
   printf("\ndata/sim (norm 30-40):\n"); for (int z = 0; z < NZ; ++z) { printf("  %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int j = 0; j < NT; ++j) { int c = z * NT + j; printf(" %6.2f", s[c] > 0 ? d[c] / (s[c] * dn / sn) : -1.); } printf("\n"); }
   for (size_t i = 0; i < iv.size(); ++i) { printf("sim median %s; zend>5 mm fraction after the table if path:\n", V[i]); for (int z = 0; z < NZ; ++z) { printf("  %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int j = 0; j < NT; ++j) { std::vector<double> x; for (auto &t : tr) if (t.cell == z * NT + j && t.x[iv[i]] > -9000) x.push_back(t.x[iv[i]]); std::sort(x.begin(), x.end()); printf(" %6.0f", x.empty() ? -1. : x[x.size() / 2]); } printf("\n"); } }
   int iz = std::find(N.begin(), N.end(), "zend") - N.begin(); printf("fraction stopping in the gas (zend > 5 mm):\n"); for (int z = 0; z < NZ; ++z) { printf("  %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int j = 0; j < NT; ++j) { double n = 0, m = 0; for (auto &t : tr) if (t.cell == z * NT + j) { ++n; if (t.x[iz] > 5) ++m; } printf(" %6.2f", n > 0 ? m / n : -1.); } printf("\n"); } }
void analyse(const char *chain, const std::vector<T> &tr0, const double *d0, const std::vector<std::string> &N, int mask, bool full) {
   std::vector<T> tr; double d[NC]; int nc = 0; for (int c = 0; c < NC; ++c) { d[c] = masked(c, mask) ? 0 : d0[c]; if (!masked(c, mask)) ++nc; } for (auto &t : tr0) if (!masked(t.cell, mask)) tr.push_back(t); printf("\n#### mask %d (%s): %d cells, ndf with X %d\n", mask, mask == 0 ? "all" : mask == 1 ? "no 10-20 column" : mask == 2 ? "no 10-20, vz 800-900 cell" : "20-35 only", nc, nc - (mask == 1 ? 4 : mask == 3 ? 3 : 5) - (KQ - 1));
   printf("\n================ %s: %zu sim tracks in cells, %.0f data ================\n", chain, tr.size(), std::accumulate(d, d + NC, 0.));
   double s0[NC][KQ] = {{0}}; for (auto &t : tr) s0[t.cell][0]++;
   Res R0 = fit(d, s0, 1); printf("angle only (band ignored): deviance %.1f\n", R0.D);
   struct O { std::string n; Res r; double med; }; std::vector<O> all;
   for (size_t v = 0; v < N.size(); ++v) { std::vector<double> xs; for (auto &t : tr) if (t.x[v] > -9000) xs.push_back(t.x[v]); if (xs.size() < tr.size() * 0.8) { printf("  (skip %s: only %zu of %zu)\n", N[v].c_str(), xs.size(), tr.size()); continue; } std::sort(xs.begin(), xs.end()); O o; o.n = N[v]; for (int q = 0; q <= KQ; ++q) o.r.q[q] = xs[std::min(xs.size() - 1, xs.size() * q / KQ)];
      double s[NC][KQ] = {{0}}; for (auto &t : tr) { if (t.x[v] <= -9000) continue; int k = 0; while (k < KQ - 1 && t.x[v] >= o.r.q[k + 1]) ++k; s[t.cell][k]++; } double q[KQ + 1]; std::copy(o.r.q, o.r.q + KQ + 1, q); o.r = fit(d, s, KQ); std::copy(q, q + KQ + 1, o.r.q); all.push_back(o); }
   std::sort(all.begin(), all.end(), [](const O &a, const O &b) { return a.r.D < b.r.D; });
   printf("%-14s %8s  g by quantile of X (edges)\n", "variable", "dev/26");
   int shown = 0; for (auto &o : all) { if (!full && ++shown > 8 && o.n != "random" && o.n != "tvz") continue; printf("%-14s %8.1f  ", o.n.c_str(), o.r.D); for (int k = 0; k < KQ; ++k) printf("%.2f ", o.r.g[k]); printf(" [%.3g", o.r.q[0]); for (int k = 1; k <= KQ; ++k) printf(" %.3g", o.r.q[k]); printf("]\n"); }
   for (int b = 0; full && b < 2 && b < (int)all.size(); ++b) { printf("\nbest %d: %s, F(theta_cm) relative to 35-40:", b + 1, all[b].n.c_str()); for (int j = 0; j < NT; ++j) printf(" %.2f", all[b].r.F[j] / all[b].r.F[NT - 1]); printf("\n data(model) by vz row, theta columns\n"); for (int z = 0; z < NZ; ++z) { printf("  %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int j = 0; j < NT; ++j) printf("  %3.0f(%5.1f)", d[z * NT + j], all[b].r.mu[z * NT + j]); printf("\n"); } }
}
}
void geom_corr_dp() { using namespace gc; gErrorIgnoreLevel = kError; gSystem->Load("libAtReconstruction.so");
   // truth variables per event
   std::vector<std::string> N = {"random", "tvz", "tth", "tke", "R0", "pitch", "nturnz", "turns", "path", "zend", "rend", "rmax", "zrmax", "zapex", "rapex", "zret", "rret", "dzturn"};
   for (int i = 0; i < 10; ++i) N.push_back(Form("zone%d", i)); const char *RN[7] = {"ring0-30", "ring30-60", "ring60-100", "ring100-150", "ring150-200", "ring200-250", "ring250+"}; for (int i = 0; i < 7; ++i) N.push_back(RN[i]);
   N.insert(N.end(), {"nh", "ncl", "qs", "zext", "rzmax", "rzmin", "d0"});
   std::map<int, std::vector<double>> X; std::map<std::pair<long, long>, int> key; TRandom3 RND(99);
   { TFile f("/mnt/f/a1975/sim_rebuild/dp_geom/truth_vars.root"); auto *t = (TTree *)f.Get("g"); int ev; double v[18], zone[10], ring[7]; t->SetBranchAddress("ev", &ev); const char *TV[17] = {"tvz", "tth", "tke", "R0", "pitch", "nturnz", "turns", "path", "zend", "rend", "rmax", "zrmax", "zapex", "rapex", "zret", "rret", "dzturn"}; for (int i = 0; i < 17; ++i) t->SetBranchAddress(TV[i], &v[i + 1]); t->SetBranchAddress("zone", zone); t->SetBranchAddress("ring", ring);
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); std::vector<double> x; x.push_back(RND.Uniform()); for (int k = 1; k < 18; ++k) x.push_back(v[k]); for (int k = 0; k < 10; ++k) x.push_back(zone[k]); for (int k = 0; k < 7; ++k) x.push_back(ring[k]); for (int k = 0; k < 7; ++k) x.push_back(-9999); X[ev] = x; key[{lround(v[1] * 1e4), lround(v[2] * 1e4)}] = ev; } }
   const size_t RECO0 = N.size() - 7;
   // sim census: reco variables, and the Spyral sim cells
   std::vector<T> TS; double dS[NC] = {0}, dG[NC] = {0};
   { TFile f("/mnt/f/a1975/sim_rebuild/dp_census/p_gs_census.root"); auto *t = (TTree *)f.Get("c"); int ev, nh, ncl, sv; double qs, zmin, zmax, rzmin, rzmax, d0, sth, svz, ske, sex, stcm; t->SetBranchAddress("ev", &ev); t->SetBranchAddress("nh", &nh); t->SetBranchAddress("ncl", &ncl); t->SetBranchAddress("sv", &sv); t->SetBranchAddress("qs", &qs); t->SetBranchAddress("zmin", &zmin); t->SetBranchAddress("zmax", &zmax); t->SetBranchAddress("rzmin", &rzmin); t->SetBranchAddress("rzmax", &rzmax); t->SetBranchAddress("d0", &d0); t->SetBranchAddress("sth", &sth); t->SetBranchAddress("svz", &svz); t->SetBranchAddress("ske", &ske); t->SetBranchAddress("sex", &sex); t->SetBranchAddress("stcm", &stcm);
     std::set<int> seen; for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); if (!sv || sth < 89 || sth > 161 || ske < 1.2 || sex < -1.5 || sex > 0.729 || !seen.insert(ev).second) continue; auto it = X.find(ev); if (it == X.end()) continue; auto &x = it->second; x[RECO0] = nh; x[RECO0 + 1] = ncl; x[RECO0 + 2] = qs; x[RECO0 + 3] = zmax - zmin; x[RECO0 + 4] = rzmax; x[RECO0 + 5] = rzmin; x[RECO0 + 6] = d0; int a = tb(stcm), z = zb(svz); if (a >= 0 && z >= 0) TS.push_back({z * NT + a, x}); } }
   // Spyral data, with the proton gate
   { TGraph G; std::ifstream in("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/proton_d2_dv1104.json"); std::string str((std::istreambuf_iterator<char>(in)), {}); auto p = str.find("vertices"); p = str.find('[', p) + 1; while (true) { auto a = str.find('[', p); if (a == std::string::npos) break; auto b = str.find(']', a); double x, y; sscanf(str.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &x, &y); G.AddPoint(x, y); p = b + 1; if (str[str.find_first_not_of(" \n\t,", p)] == ']') break; }
     TChain c("c"); c.Add("/mnt/f/a1975/census_arc/data_*.root"); int sv; double sth, svz, ske, sex, stcm, ssq, sbr; c.SetBranchAddress("sv", &sv); c.SetBranchAddress("sth", &sth); c.SetBranchAddress("svz", &svz); c.SetBranchAddress("ske", &ske); c.SetBranchAddress("sex", &sex); c.SetBranchAddress("stcm", &stcm); c.SetBranchAddress("ssq", &ssq); c.SetBranchAddress("sbr", &sbr); c.SetBranchStatus("*", 0); for (auto b : {"sv", "sth", "svz", "ske", "sex", "stcm", "ssq", "sbr"}) c.SetBranchStatus(b, 1);
     for (Long64_t i = 0; i < c.GetEntries(); ++i) { c.GetEntry(i); if (!sv || sth < 89 || sth > 161 || ske < 1.2 || sex < -1.5 || sex > 0.729 || !G.IsInside(ssq, sbr)) continue; int a = tb(stcm), z = zb(svz); if (a >= 0 && z >= 0) dS[z * NT + a]++; } }
   // genfit, frame B
   std::vector<T> TG;
   { TFile f("/mnt/f/a1975/sim_rebuild/dp_eps/new_dedup.root"); auto *t = (TTree *)f.Get("r"); double ke, th, vz, c2, tth, tvz; t->SetBranchAddress("ke", &ke); t->SetBranchAddress("th", &th); t->SetBranchAddress("vz", &vz); t->SetBranchAddress("c2", &c2); t->SetBranchAddress("tth", &tth); t->SetBranchAddress("tvz", &tvz); int miss = 0;
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); double ex, tc; if (!(c2 <= 13.5) || th < 89 || th > 161 || ke < 1.2 || !kine2b(EbB(vz), th * TMath::DegToRad(), ke, ex, tc) || ex < -1.5 || ex > 0.729) continue; int a = tb(tc), z = zb(vz); if (a < 0 || z < 0) continue; auto k = key.find({lround(tvz * 1e4), lround(tth * 1e4)}); if (k == key.end()) { ++miss; continue; } TG.push_back({z * NT + a, X[k->second]}); } printf("genfit sim tracks without truth match: %d\n", miss); }
   { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); auto *t = (TTree *)f.Get("pk"); float ke, th, vz, c2; t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vertexz", &vz); t->SetBranchAddress("chi2ndf", &c2); for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); double ex, tc; if (!(c2 <= 13.5) || th < 89 || th > 161 || ke < 1.2 || !kine2b(EbB(vz), th * TMath::DegToRad(), ke, ex, tc) || ex < -1.5 || ex > 0.729) continue; int a = tb(tc), z = zb(vz); if (a >= 0 && z >= 0) dG[z * NT + a]++; } }
   for (int m = 3; m <= 3; ++m) { analyse("GENFIT frame B", TG, dG, N, m, false); analyse("SPYRAL census", TS, dS, N, m, false); }
}
