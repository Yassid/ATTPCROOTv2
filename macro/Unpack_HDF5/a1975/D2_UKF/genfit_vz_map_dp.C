#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
/// @file genfit_vz_map_dp.C  The loss-band table from the adopted genfit chain: data / sim in (vz, theta_cm), counts in the bound window, sim fit at -2.85 (new_dedup). Cuts of angdist_bound_dp.C; frame B and frame A beam energy.
namespace gm { using namespace exsel; const int NT = 5, NZ = 8; const double TB[NT + 1] = {10, 20, 25, 30, 35, 40}, ZB[NZ + 1] = {100, 200, 300, 400, 500, 600, 700, 800, 900};
double Eb(int fr, double vz) { return fr ? 183.873 - 0.011099 * vz : 184.25 - 0.011030 * (vz - 500); }
int tb(double t) { for (int i = 0; i < NT; ++i) if (t >= TB[i] && t < TB[i + 1]) return i; return -1; } int zb(double z) { for (int i = 0; i < NZ; ++i) if (z >= ZB[i] && z < ZB[i + 1]) return i; return -1; } }
void genfit_vz_map_dp() { using namespace gm;
   struct R { double ke, th, vz, c2; }; std::vector<R> sim, dat;
   { TFile f("/mnt/f/a1975/sim_rebuild/dp_eps/new_dedup.root"); auto *t = (TTree *)f.Get("r"); R s; t->SetBranchAddress("ke", &s.ke); t->SetBranchAddress("th", &s.th); t->SetBranchAddress("vz", &s.vz); t->SetBranchAddress("c2", &s.c2); for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); sim.push_back(s); } }
   { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); auto *t = (TTree *)f.Get("pk"); float ke, th, vz, c2; t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vertexz", &vz); t->SetBranchAddress("chi2ndf", &c2); for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); dat.push_back({ke, th, vz, c2}); } }
   for (int fr = 0; fr < 2; ++fr) { double Nd[NZ][NT] = {{0}}, Ns[NZ][NT] = {{0}}; const double TOP = 0.729;
      auto fill = [&](const std::vector<R> &v, double (*N)[NT]) { for (auto &d : v) { double ex, tc; if (!(d.c2 <= 13.5) || d.th < 89 || d.th > 161 || d.ke < 1.2 || !kine2b(Eb(fr, d.vz), d.th * TMath::DegToRad(), d.ke, ex, tc) || ex < -1.5 || ex > TOP) continue; int z = zb(d.vz), t = tb(tc); if (z >= 0 && t >= 0) N[z][t]++; } };
      fill(dat, Nd); fill(sim, Ns); double dn = 0, sn = 0; for (int z = 1; z < 8; ++z) for (int t = 3; t < NT; ++t) { dn += Nd[z][t]; sn += Ns[z][t]; }
      printf("\n=== frame %c, norm data/sim on 30-40, vz 200-900 = %.4f; data counts (sim counts); ratio +- stat ===\n vz        ", fr ? 'B' : 'A', dn / sn); for (int t = 0; t < NT; ++t) printf("        %2.0f-%2.0f      ", TB[t], TB[t + 1]); printf("\n");
      for (int z = 0; z < NZ; ++z) { printf(" %3.0f-%3.0f ", ZB[z], ZB[z + 1]); for (int t = 0; t < NT; ++t) { double d = Nd[z][t], s = Ns[z][t]; if (s > 0) printf(" %4.0f(%3.0f) %4.2f+-%4.2f", d, s, d / (s * dn / sn), sqrt(std::max(d, 1.)) / (s * dn / sn)); else printf(" %4.0f(  0)     -     ", d); } printf("\n"); } } }
