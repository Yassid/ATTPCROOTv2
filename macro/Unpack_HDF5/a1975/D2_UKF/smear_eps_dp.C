#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
/// @file smear_eps_dp.C
/// (d,p) bound-state efficiency with the data's Ex width added to the sim, at the right field sign. Reads the dump of sim_dump_dp.C.
/// Same selection and window as angdist_bound_dp.C; vz 200-400; truth window on z_true + 19.1 mm.
///   root -l -b -q 'smear_eps_dp.C("dump.root","outdir")'
namespace sm {
using namespace exsel;
const int NT = 7; const double TB[NT + 1] = {10, 15, 20, 25, 30, 35, 40, 45};
const double EXLO = -1.5, EXHI = 0.729, Z0 = 200, Z1 = 400, LW = 316.98 * 200 / 971.73, ZOFF = 19.1;
const int NTOY = 400;
double Eb(double vz) { return 184.25 - 0.011030 * (vz - 500); }
bool sel(double ke, double th, double vz, double c2, double &ex, double &tc) {
   if (!(c2 <= 13.5) || !(th >= 89 && th <= 161) || !(ke >= 1.2)) return false;
   return kine2b(Eb(vz), th * TMath::DegToRad(), ke, ex, tc); }
int tb(double t) { for (int i = 0; i < NT; ++i) if (t >= TB[i] && t < TB[i + 1]) return i; return -1; }
struct Fit { double n, mu, emu, sg, esg, bg, ebg; int st; };
/// gaussian + flat, Poisson likelihood, inside -4..window top only
Fit fitpeak(TH1D *h) {
   TF1 f("f", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -4, 0.7);
   f.SetParameters(h->GetMaximum(), -0.2, 0.3, 0.5); f.SetParLimits(0, 0, 1e6); f.SetParLimits(1, -1.5, 0.7); f.SetParLimits(2, 0.03, 1.5); f.SetParLimits(3, 0, 1e5);
   int st = h->Fit(&f, "LQN0SR") ;
   return {f.GetParameter(0), f.GetParameter(1), f.GetParError(1), f.GetParameter(2), f.GetParError(2), f.GetParameter(3), f.GetParError(3), st}; }
void pf(const char *l, const Fit &f) { printf("  %-22s amp %7.1f  mu %+.3f(%.3f)  sigma %.3f(%.3f)  flat %6.2f(%.2f) %s\n", l, f.n, f.mu, f.emu, f.sg, f.esg, f.bg, f.ebg, (f.sg < 0.035 || f.sg > 1.45 || f.mu < -1.45 || f.mu > 0.65 || f.bg < 1e-3) ? "RAILED?" : ""); }
struct S { double ke, th, vz, c2, tke, tth, tvz; };
struct D { double ke, th, vz, c2; };
}
void smear_eps_dp(TString dump = "/mnt/f/a1975/sim_rebuild/dp_eps/new_dedup.root", TString outdir = "/mnt/f/a1975/sim_rebuild/dp_smear") {
   using namespace sm;
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   gErrorIgnoreLevel = kError; ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
   std::vector<std::array<double, 3>> gen;
   { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root"); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *mc = nullptr; t->SetBranchStatus("*", 0); t->SetBranchStatus("MCTrack*", 1); t->SetBranchAddress("MCTrack", &mc);
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i);
        for (int j = 0; j < mc->GetEntriesFast(); ++j) { auto *m = (AtMCTrack *)mc->At(j); if (!m || m->GetPdgCode() != 2212) continue;
           double p = m->GetP() * 1000; gen.push_back({m->GetStartZ() * 10.0, acos(m->GetPz() / m->GetP()) * TMath::RadToDeg(), sqrt(p * p + 938.272 * 938.272) - 938.272}); break; } } }
   std::vector<S> sim;
   { TFile f(dump); auto *t = (TTree *)f.Get("r"); S s; t->SetBranchAddress("ke", &s.ke); t->SetBranchAddress("th", &s.th); t->SetBranchAddress("vz", &s.vz); t->SetBranchAddress("c2", &s.c2); t->SetBranchAddress("tke", &s.tke); t->SetBranchAddress("tth", &s.tth); t->SetBranchAddress("tvz", &s.tvz);
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); sim.push_back(s); } }
   std::vector<D> dat;
   { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); auto *t = (TTree *)f.Get("pk"); float ke, th, vz, c2; t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vertexz", &vz); t->SetBranchAddress("chi2ndf", &c2);
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); dat.push_back({ke, th, vz, c2}); } }
   printf("generated %zu  sim fitted %zu  data %zu\n", gen.size(), sim.size(), dat.size());
   // Ng: truth window in the frame of the reco vz
   double Ng[NT] = {0}, Y[NT] = {0}, Bg[NT] = {0};
   for (auto &g : gen) { double z = g[0] + ZOFF; if (z < Z0 || z >= Z1) continue; double ex, tc; if (!kine2b(Eb(z), g[1] * TMath::DegToRad(), g[2], ex, tc)) continue; int b = tb(tc); if (b >= 0) Ng[b]++; }
   // peaks: 4 bins + pooled, data at vz 200-400 and 200-900, sim at 200-400
   const int NP = 5; const double PL[NP] = {20, 25, 30, 35, 20}, PH[NP] = {25, 30, 35, 40, 40};
   TH1D *hd[NP], *hw[NP], *hs[NP];
   for (int p = 0; p < NP; ++p) { hd[p] = new TH1D(Form("hd%d", p), "", 80, -5, 3); hw[p] = new TH1D(Form("hw%d", p), "", 80, -5, 3); hs[p] = new TH1D(Form("hs%d", p), "", 80, -5, 3); }
   for (auto &d : dat) { double ex, tc; if (d.vz < Z0 || d.vz >= 900) continue; if (!sel(d.ke, d.th, d.vz, d.c2, ex, tc)) continue;
      for (int p = 0; p < NP; ++p) if (tc >= PL[p] && tc < PH[p]) { hw[p]->Fill(ex); if (d.vz < Z1) hd[p]->Fill(ex); }
      if (d.vz < Z1) { int b = tb(tc); if (b >= 0) { if (ex >= EXLO && ex <= EXHI) Y[b]++; else if (ex >= -4 && ex < EXLO) Bg[b]++; } } }
   std::vector<S> ss; std::vector<double> sex, stc;
   for (auto &s : sim) { double ex, tc; if (s.vz < Z0 || s.vz >= Z1) continue; if (!sel(s.ke, s.th, s.vz, s.c2, ex, tc)) continue; ss.push_back(s); sex.push_back(ex); stc.push_back(tc);
      for (int p = 0; p < NP; ++p) if (tc >= PL[p] && tc < PH[p]) hs[p]->Fill(ex); }
   Fit fd[NP], fw[NP], fs[NP];
   printf("\n--- Ex peak fits (gaus + flat, -4..0.7, Poisson L, 100 keV) ---\n");
   for (int p = 0; p < NP; ++p) { fd[p] = fitpeak(hd[p]); fw[p] = fitpeak(hw[p]); fs[p] = fitpeak(hs[p]); printf(" theta_cm %.0f-%.0f\n", PL[p], PH[p]); pf("data vz 200-400", fd[p]); pf("data vz 200-900", fw[p]); pf("sim  vz 200-400", fs[p]); }
   // response of the reco theta_cm
   printf("\n--- theta_cm reco - truth, sim, no Ex window ---\n");
   { std::vector<double> dv[NT]; double mig[NT][NT] = {{0}};
     for (size_t i = 0; i < ss.size(); ++i) { double ex, tc, z = ss[i].tvz + ZOFF; if (!kine2b(Eb(z), ss[i].tth * TMath::DegToRad(), ss[i].tke, ex, tc)) continue; int bt = tb(tc), br = tb(stc[i]); if (bt >= 0) dv[bt].push_back(stc[i] - tc); if (bt >= 0 && br >= 0) mig[bt][br]++; }
     for (int b = 0; b < NT; ++b) { if (dv[b].size() < 5) continue; std::sort(dv[b].begin(), dv[b].end()); size_t n = dv[b].size(); printf("  true %2.0f-%2.0f  n %4zu  p16 %+.2f  med %+.2f  p84 %+.2f   reco-bin fractions:", TB[b], TB[b + 1], n, dv[b][n * 16 / 100], dv[b][n / 2], dv[b][n * 84 / 100]); double sum = 0; for (int r = 0; r < NT; ++r) sum += mig[b][r]; for (int r = 0; r < NT; ++r) printf(" %.3f", sum > 0 ? mig[b][r] / sum : 0.); printf("\n"); } }
   // variants: added smear from the pooled peak (A), per bin (B); shift to data centroid (C); then same width carried by theta_lab or KE, which also moves theta_cm
   const double s0 = sqrt(std::max(0., fd[4].sg * fd[4].sg - fs[4].sg * fs[4].sg));
   double sb[NT] = {0}, sh[NT] = {0}; for (int b = 0; b < NT; ++b) { sb[b] = s0; for (int p = 0; p < 4; ++p) if (TB[b] == PL[p]) { sb[b] = sqrt(std::max(0., fd[p].sg * fd[p].sg - fs[p].sg * fs[p].sg)); sh[b] = fd[p].mu - fs[p].mu; } }
   printf("\n--- added smear (MeV): pooled %.3f; per bin", s0); for (int b = 2; b < 6; ++b) printf(" %.3f", sb[b]); printf("; shift per bin"); for (int b = 2; b < 6; ++b) printf(" %+.3f", sh[b]); printf("; pooled shift %+.3f\n", fd[4].mu - fs[4].mu);
   const int NV = 7; const char *VN[NV] = {"none", "Ex_pool", "Ex_perbin", "Ex_pool_shift", "Ex_perbin_shift", "thlab_pool", "ke_pool"};
   double Nr[NV][NT] = {{0}}; TRandom3 R(4242);
   for (size_t i = 0; i < ss.size(); ++i) { const S &s = ss[i]; int b0 = tb(stc[i]);
      if (b0 >= 0 && sex[i] >= EXLO && sex[i] <= EXHI) Nr[0][b0] += 1;
      if (b0 >= 0) for (int t = 0; t < NTOY; ++t) { double w = 1. / NTOY, g = R.Gaus(), e;
         e = sex[i] + s0 * g; if (e >= EXLO && e <= EXHI) Nr[1][b0] += w;
         e = sex[i] + sb[b0] * g; if (e >= EXLO && e <= EXHI) Nr[2][b0] += w;
         e = sex[i] + (fd[4].mu - fs[4].mu) + s0 * g; if (e >= EXLO && e <= EXHI) Nr[3][b0] += w;
         e = sex[i] + sh[b0] + sb[b0] * g; if (e >= EXLO && e <= EXHI) Nr[4][b0] += w; }
      double ea, tca, da = 0, dk = 0;
      if (kine2b(Eb(s.vz), (s.th + 0.1) * TMath::DegToRad(), s.ke, ea, tca)) da = (ea - sex[i]) / 0.1;
      if (kine2b(Eb(s.vz), s.th * TMath::DegToRad(), s.ke + 0.02, ea, tca)) dk = (ea - sex[i]) / 0.02;
      for (int t = 0; t < NTOY; ++t) { double w = 1. / NTOY, g = R.Gaus(), ex, tc; int b;
         if (da != 0) { double th = s.th + s0 / fabs(da) * g; if (th >= 89 && th <= 161 && kine2b(Eb(s.vz), th * TMath::DegToRad(), s.ke, ex, tc) && ex >= EXLO && ex <= EXHI && (b = tb(tc)) >= 0) Nr[5][b] += w; }
         if (dk != 0) { double ke = s.ke + s0 / fabs(dk) * g; if (ke >= 1.2 && kine2b(Eb(s.vz), s.th * TMath::DegToRad(), ke, ex, tc) && ex >= EXLO && ex <= EXHI && (b = tb(tc)) >= 0) Nr[6][b] += w; } } }
   printf("\n--- efficiency, Ng"); for (int b = 0; b < NT; ++b) printf(" %.0f", Ng[b]); printf(" ---\n%-18s", "theta_cm"); for (int b = 0; b < NT; ++b) printf("  %2.0f-%2.0f", TB[b], TB[b + 1]); printf("\n");
   for (int v = 0; v < NV; ++v) { printf("%-18s", VN[v]); for (int b = 0; b < NT; ++b) printf("  %.3f", Ng[b] > 0 ? Nr[v][b] / Ng[b] : 0.); printf("\n"); }
   for (int v = 0; v < NV; ++v) { FILE *fo = fopen(outdir + Form("/ang_%s.csv", VN[v]), "w"); fprintf(fo, "# smear_eps_dp.C %s, variant %s; vz 200-400, truth window z+19.1; NOT ADOPTED\nthcm_lo,thcm_hi,thcm_mid,Y,bg,eps,eps_err,dsdo_mb_sr,err_mb_sr,use\n", dump.Data(), VN[v]);
      printf("dsdo %-13s", VN[v]);
      for (int b = 0; b < NT; ++b) { double e = Ng[b] > 0 ? Nr[v][b] / Ng[b] : 0, ee = Ng[b] > 0 ? sqrt(std::max(e * fabs(1 - e), 1. / Ng[b] / Ng[b]) / Ng[b]) : 0;
         double dO = 2 * TMath::Pi() * (cos(TB[b] * TMath::DegToRad()) - cos(TB[b + 1] * TMath::DegToRad())), bgs = Bg[b] * (2.229 / 2.5), net = Y[b] - bgs, x = e > 0 ? net / (e * dO * LW) : 0;
         double rel2 = net > 0 ? (Y[b] + bgs * bgs / std::max(1., Bg[b])) / (net * net) : 0; rel2 += e > 0 ? (ee / e) * (ee / e) : 0; double er = x * sqrt(rel2); int use = b >= 2 && b <= 5;
         fprintf(fo, "%.0f,%.0f,%.1f,%.0f,%.2f,%.4f,%.4f,%.4f,%.4f,%d\n", TB[b], TB[b + 1], 0.5 * (TB[b] + TB[b + 1]), Y[b], bgs, e, ee, x, er, use); printf(" %5.2f(%.2f)", x, er); }
      printf("\n"); fclose(fo); }
   // plots of the sim and data Ex peaks
   gStyle->SetOptStat(0); TCanvas c("c", "", 1500, 800); c.Divide(3, 2);
   for (int p = 0; p < NP; ++p) { c.cd(p + 1); gPad->SetGrid(); hd[p]->SetTitle(Form("#theta_{cm} %.0f-%.0f, vz 200-400;E_{x} [MeV];counts / 100 keV", PL[p], PH[p])); hd[p]->GetXaxis()->SetRangeUser(-5, 3); hd[p]->SetLineColor(kBlue + 1); hd[p]->SetLineWidth(2);
      double sc = hs[p]->Integral(hs[p]->FindBin(-1.5), hs[p]->FindBin(0.7)) > 0 ? hd[p]->Integral(hd[p]->FindBin(-1.5), hd[p]->FindBin(0.7)) / hs[p]->Integral(hs[p]->FindBin(-1.5), hs[p]->FindBin(0.7)) : 1; hs[p]->Scale(sc); hs[p]->SetLineColor(kRed + 1); hd[p]->SetMaximum(1.15 * std::max(hd[p]->GetMaximum(), hs[p]->GetMaximum())); hd[p]->Draw("E"); hs[p]->Draw("hist same"); }
   c.SaveAs(outdir + "/peaks.png");
}
