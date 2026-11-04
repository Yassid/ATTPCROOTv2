#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
/// @file window_check_dp.C
/// (d,p) bound window: loss per bound state and per angle bin, in the frame of angdist_bound_dp.C (A, pivot Ebeam) and in the frame where a generated g.s. reads Ex = 0 (B, Ebeam from sim truth).
/// Adds Ex to the g.s. sim to stand in for 0.217 and 0.331. Also Y vs window top, and a level-mix fit with the smear free.
///   root -l -b -q 'window_check_dp.C("dump.root","outdir")'
namespace wc {
using namespace exsel;
const int NT = 7; const double TB[NT + 1] = {10, 15, 20, 25, 30, 35, 40, 45};
const double Z0 = 200, Z1 = 400, LW = 316.98 * 200 / 971.73, ZOFF = 19.1, EXLO = -1.5;
double gA = 0, gB = 0;
double EbA(double vz) { return 184.25 - 0.011030 * (vz - 500); }
double EbB(double vz) { return gA + gB * vz; }
double Eb(int fr, double vz) { return fr == 0 ? EbA(vz) : EbB(vz); }
int tb(double t) { for (int i = 0; i < NT; ++i) if (t >= TB[i] && t < TB[i + 1]) return i; return -1; }
bool sel(int fr, double ke, double th, double vz, double c2, double &ex, double &tc) {
   if (!(c2 <= 13.5) || !(th >= 89 && th <= 161) || !(ke >= 1.2)) return false;
   return kine2b(Eb(fr, vz), th * TMath::DegToRad(), ke, ex, tc); }
struct Fit { double n, mu, emu, sg, esg, bg; };
Fit fitpeak(TH1D *h, double top) {
   TF1 f("f", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -4, top);
   f.SetParameters(h->GetMaximum(), -0.1, 0.3, 0.5); f.SetParLimits(0, 0, 1e6); f.SetParLimits(1, -1.5, 0.7); f.SetParLimits(2, 0.03, 1.5); f.SetParLimits(3, 0, 1e5);
   h->Fit(&f, "LQN0SR");
   return {f.GetParameter(0), f.GetParameter(1), f.GetParError(1), f.GetParameter(2), f.GetParError(2), f.GetParameter(3)}; }
struct S { double ke, th, vz, c2, tke, tth, tvz; };
struct D { double ke, th, vz, c2; };
/// bin fractions of the sim line shifted by d and smeared by s, on -1.5..1.0, 50 keV
const int NB = 50; const double XL = -1.5, XW = 0.05;
void tmpl(TH1D *fine, double d, double s, double *t) {
   for (int j = 0; j < NB; ++j) t[j] = 0;
   for (int i = 1; i <= fine->GetNbinsX(); ++i) { double c = fine->GetBinContent(i); if (c <= 0) continue; double m = fine->GetBinCenter(i) + d;
      for (int j = 0; j < NB; ++j) { double a = XL + j * XW; t[j] += c * 0.5 * (erf((a + XW - m) / (s * 1.41421356)) - erf((a - m) / (s * 1.41421356))); } }
   double tot = fine->Integral(); for (int j = 0; j < NB; ++j) t[j] = t[j] / tot + 1e-12; }
/// EM on a Poisson mixture: templates fixed, amplitudes >= 0. Returns -2lnL
double em(TH1D *h, int nt, double T[][NB], double *a) {
   double n[NB], tot = 0; for (int j = 0; j < NB; ++j) { n[j] = h->GetBinContent(h->FindBin(XL + (j + 0.5) * XW)); tot += n[j]; }
   for (int k = 0; k < nt; ++k) a[k] = tot / nt;
   for (int it = 0; it < 300; ++it) { double na[4] = {0}; for (int j = 0; j < NB; ++j) { double m = 0; for (int k = 0; k < nt; ++k) m += a[k] * T[k][j]; if (n[j] > 0) for (int k = 0; k < nt; ++k) na[k] += n[j] * a[k] * T[k][j] / m; } for (int k = 0; k < nt; ++k) a[k] = na[k]; }
   double L = 0; for (int j = 0; j < NB; ++j) { double m = 0; for (int k = 0; k < nt; ++k) m += a[k] * T[k][j]; L += 2 * (m - n[j] * log(m)); if (n[j] > 0) L += 2 * (n[j] * log(n[j]) - n[j]); }
   return L; }
void scan(const char *lab, TH1D *h, TH1D *fine, bool three, bool freeD) {
   double best = 1e30, bs = 0, bd = 0, ba[4] = {0}; static double T[4][NB]; double a[4]; const double E[3] = {0, 0.217, 0.331}; int nl = three ? 3 : 1;
   for (double s = 0.03; s <= 0.45; s += 0.02) for (double d = freeD ? -0.34 : 0; d <= (freeD ? 0.2001 : 0); d += 0.02) {
      for (int k = 0; k < nl; ++k) tmpl(fine, d + E[k], s, T[k]);
      for (int j = 0; j < NB; ++j) T[nl][j] = 1. / NB;
      double L = em(h, nl + 1, T, a); if (L < best) { best = L; bs = s; bd = d; for (int k = 0; k <= nl; ++k) ba[k] = a[k]; } }
   double sum = 0; for (int k = 0; k < nl; ++k) sum += ba[k];
   printf("  %-8s %-24s -2lnL %7.2f  smear %.2f  shift %+.2f  N %6.1f flat %5.1f", lab, three ? (freeD ? "3 levels, shift free" : "3 levels, shift 0") : (freeD ? "g.s., shift free" : "g.s., shift 0"), best, bs, bd, sum, ba[nl]);
   if (three) printf("  frac %.2f %.2f %.2f", ba[0] / sum, ba[1] / sum, ba[2] / sum); printf("\n"); }
}
void window_check_dp(TString dump = "/mnt/f/a1975/sim_rebuild/dp_eps/new_dedup.root", TString outdir = "/mnt/f/a1975/sim_rebuild/dp_window") {
   using namespace wc;
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so"); gErrorIgnoreLevel = kError; ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
   std::vector<std::array<double, 3>> gen;
   { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root"); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *mc = nullptr; t->SetBranchStatus("*", 0); t->SetBranchStatus("MCTrack*", 1); t->SetBranchAddress("MCTrack", &mc);
     for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); for (int j = 0; j < mc->GetEntriesFast(); ++j) { auto *m = (AtMCTrack *)mc->At(j); if (!m || m->GetPdgCode() != 2212) continue; double p = m->GetP() * 1000; gen.push_back({m->GetStartZ() * 10.0, acos(m->GetPz() / m->GetP()) * TMath::RadToDeg(), sqrt(p * p + 938.272 * 938.272) - 938.272}); break; } } }
   std::vector<S> sim; { TFile f(dump); auto *t = (TTree *)f.Get("r"); S s; t->SetBranchAddress("ke", &s.ke); t->SetBranchAddress("th", &s.th); t->SetBranchAddress("vz", &s.vz); t->SetBranchAddress("c2", &s.c2); t->SetBranchAddress("tke", &s.tke); t->SetBranchAddress("tth", &s.tth); t->SetBranchAddress("tvz", &s.tvz); for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); sim.push_back(s); } }
   std::vector<D> dat; { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); auto *t = (TTree *)f.Get("pk"); float ke, th, vz, c2; t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vertexz", &vz); t->SetBranchAddress("chi2ndf", &c2); for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); dat.push_back({ke, th, vz, c2}); } }
   // frame B: beam energy zeroing the Ex of the truth, per event, against the reco frame vz = z_true + 19.1
   { TGraph g; int nn = 0; for (auto &q : gen) { if (q[1] < 89 || q[1] > 161 || (++nn) % 3) continue; double lo = 170, hi = 200, ex, tc; bool ok = true; for (int it = 0; it < 30 && ok; ++it) { double m = 0.5 * (lo + hi); ok = kine2b(m, q[1] * TMath::DegToRad(), q[2], ex, tc); if (ex > 0) hi = m; else lo = m; } if (ok) g.AddPoint(q[0] + ZOFF, 0.5 * (lo + hi)); }
     TF1 l("l", "pol1", 100, 900); g.Fit(&l, "QNRS"); gA = l.GetParameter(0); gB = l.GetParameter(1); printf("frame B: Ebeam(vz) = %.3f %+.6f vz  (from %d truth points); at vz 300 frame A %.2f, frame B %.2f\n", gA, gB, g.GetN(), EbA(300), EbB(300)); }
   // sensitivities in the data frame
   printf("\ndEx per unit, frame A, sim events (median): theta_cm  dEx/dthlab [MeV/deg]  dEx/dKE [MeV/MeV]  dEx/dEbeam\n");
   { std::vector<double> a[NT], k[NT], e[NT]; for (auto &s : sim) { double ex, tc, x, t2; if (s.vz < Z0 || s.vz >= Z1 || !sel(0, s.ke, s.th, s.vz, s.c2, ex, tc)) continue; int b = tb(tc); if (b < 0) continue; if (kine2b(EbA(s.vz), (s.th + 0.1) * TMath::DegToRad(), s.ke, x, t2)) a[b].push_back((x - ex) / 0.1); if (kine2b(EbA(s.vz), s.th * TMath::DegToRad(), s.ke + 0.02, x, t2)) k[b].push_back((x - ex) / 0.02); if (kine2b(EbA(s.vz) + 1, s.th * TMath::DegToRad(), s.ke, x, t2)) e[b].push_back(x - ex); }
     for (int b = 2; b < 6; ++b) { std::sort(a[b].begin(), a[b].end()); std::sort(k[b].begin(), k[b].end()); std::sort(e[b].begin(), e[b].end()); printf("   %2.0f-%2.0f   %+.3f  %+.3f  %+.4f\n", TB[b], TB[b + 1], a[b][a[b].size() / 2], k[b][k[b].size() / 2], e[b][e[b].size() / 2]); } }
   const int NP = 5; const double PL[NP] = {20, 25, 30, 35, 20}, PH[NP] = {25, 30, 35, 40, 40};
   const double TOP[3] = {0.729, 0.5, 1.0}, LEV[3] = {0, 0.217, 0.331};
   for (int fr = 0; fr < 2; ++fr) {
      printf("\n================ FRAME %s ================\n", fr ? "B (generated g.s. at Ex 0)" : "A (angdist_bound_dp.C, pivot)");
      TH1D *hd[NP], *hs[NP], *hf[NP]; double Ng[NT] = {0}, Y[3][NT] = {{0}}, Bg[NT] = {0};
      for (int p = 0; p < NP; ++p) { hd[p] = new TH1D(Form("hd%d_%d", fr, p), "", 160, -5, 3); hs[p] = new TH1D(Form("hs%d_%d", fr, p), "", 160, -5, 3); hf[p] = new TH1D(Form("hf%d_%d", fr, p), "", 500, -2.5, 2.5); }
      for (auto &g : gen) { double z = g[0] + ZOFF, ex, tc; if (z < Z0 || z >= Z1 || !kine2b(Eb(fr, z), g[1] * TMath::DegToRad(), g[2], ex, tc)) continue; int b = tb(tc); if (b >= 0) Ng[b]++; }
      for (auto &d : dat) { double ex, tc; if (d.vz < Z0 || d.vz >= Z1 || !sel(fr, d.ke, d.th, d.vz, d.c2, ex, tc)) continue; for (int p = 0; p < NP; ++p) if (tc >= PL[p] && tc < PH[p]) hd[p]->Fill(ex); int b = tb(tc); if (b < 0) continue; for (int w = 0; w < 3; ++w) if (ex >= EXLO && ex <= TOP[w]) Y[w][b]++; if (ex >= -4 && ex < EXLO) Bg[b]++; }
      std::vector<double> sex; std::vector<int> sb; for (auto &s : sim) { double ex, tc; if (s.vz < Z0 || s.vz >= Z1 || !sel(fr, s.ke, s.th, s.vz, s.c2, ex, tc)) continue; int b = tb(tc); if (b < 0) continue; sex.push_back(ex); sb.push_back(b); for (int p = 0; p < NP; ++p) if (tc >= PL[p] && tc < PH[p]) { hs[p]->Fill(ex); hf[p]->Fill(ex); } }
      double top = fr ? 0.45 : 0.7; Fit fd[NP], fs[NP];
      for (int p = 0; p < NP; ++p) { fd[p] = fitpeak(hd[p], top); fs[p] = fitpeak(hs[p], top); printf(" %2.0f-%2.0f  data mu %+.3f(%.3f) sigma %.3f(%.3f) flat %5.2f | sim mu %+.3f(%.3f) sigma %.3f(%.3f) flat %5.2f\n", PL[p], PH[p], fd[p].mu, fd[p].emu, fd[p].sg, fd[p].esg, fd[p].bg, fs[p].mu, fs[p].emu, fs[p].sg, fs[p].esg, fs[p].bg); }
      const double s0 = sqrt(std::max(0., fd[4].sg * fd[4].sg - fs[4].sg * fs[4].sg)); printf(" added smear (pooled) %.3f\n", s0);
      printf(" Y by window top %.3f/%.3f/%.3f per bin 10-15..40-45:\n", TOP[0], TOP[1], TOP[2]); for (int w = 0; w < 3; ++w) { printf("   top %.3f:", TOP[w]); for (int b = 0; b < NT; ++b) printf(" %4.0f", Y[w][b]); printf("\n"); } printf("   side -4..-1.5:"); for (int b = 0; b < NT; ++b) printf(" %4.0f", Bg[b]); printf("\n");
      double ee[3][3][NT] = {{{0}}}; // [level][top][bin]
      TRandom3 R(4242 + fr); const int NTOY = 300;
      for (size_t i = 0; i < sex.size(); ++i) for (int t = 0; t < NTOY; ++t) { double g = R.Gaus(); for (int l = 0; l < 3; ++l) { double e = sex[i] + LEV[l] + s0 * g; for (int w = 0; w < 3; ++w) if (e >= EXLO && e <= TOP[w]) ee[l][w][sb[i]] += 1. / NTOY; } }
      printf(" efficiency with smear (level; bins 20-25 25-30 30-35 35-40):\n"); for (int w = 0; w < 3; ++w) for (int l = 0; l < 3; ++l) { printf("   top %.3f  Ex* %.3f:", TOP[w], LEV[l]); for (int b = 2; b < 6; ++b) printf(" %.3f", ee[l][w][b] / Ng[b]); printf("\n"); }
      for (int w = 0; w < 3; ++w) { FILE *fo = fopen(outdir + Form("/ang_fr%c_top%03.0f.csv", fr ? 'B' : 'A', TOP[w] * 1000), "w"); fprintf(fo, "# window_check_dp.C frame %c, top %.3f, eps g.s. + smear %.3f; vz 200-400; NOT ADOPTED\nthcm_lo,thcm_hi,thcm_mid,Y,bg,eps,eps_err,dsdo_mb_sr,err_mb_sr,use\n", fr ? 'B' : 'A', TOP[w], s0);
         printf(" dsdo top %.3f (g.s. eps):", TOP[w]);
         for (int b = 0; b < NT; ++b) { double e = Ng[b] > 0 ? ee[0][w][b] / Ng[b] : 0, er = Ng[b] > 0 ? sqrt(std::max(e * fabs(1 - e), 1e-6) / Ng[b]) : 0, dO = 2 * TMath::Pi() * (cos(TB[b] * TMath::DegToRad()) - cos(TB[b + 1] * TMath::DegToRad())), bgs = Bg[b] * ((TOP[w] - EXLO) / 2.5), net = Y[w][b] - bgs, x = e > 0 ? net / (e * dO * LW) : 0, rel2 = net > 0 ? (Y[w][b] + bgs * bgs / std::max(1., Bg[b])) / (net * net) : 0; rel2 += e > 0 ? er * er / e / e : 0; fprintf(fo, "%.0f,%.0f,%.1f,%.0f,%.2f,%.4f,%.4f,%.4f,%.4f,%d\n", TB[b], TB[b + 1], 0.5 * (TB[b] + TB[b + 1]), Y[w][b], bgs, e, er, x, x * sqrt(rel2), b >= 2 && b <= 5); if (b >= 2 && b <= 5) printf(" %5.2f(%.2f)", x, x * sqrt(rel2)); } printf("\n"); fclose(fo); }
      if (fr == 1) { printf("\n level-mix fits, -1.5..1.0, 50 keV, flat background, smear free\n");
         for (int p = 0; p < NP; ++p) { TString lab = Form("%.0f-%.0f", PL[p], PH[p]); hd[p]->Rebin(1); scan(lab, hd[p], hf[p], false, false); scan(lab, hd[p], hf[p], true, false); scan(lab, hd[p], hf[p], false, true); scan(lab, hd[p], hf[p], true, true); } }
   }
}
