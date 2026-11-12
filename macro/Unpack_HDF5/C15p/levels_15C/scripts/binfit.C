#pragma once
#include <Math/SpecFuncMathMore.h>
#include "binacc.h"
#include "psmodel.h"
#include "binbuild.C"
#include "lineshape.C"
// =============================================================================================
//  PER-theta_CM-BIN fit of 15C(d,d'), mirroring C15_pp_ana_AngBins_Barrier_6BW.C.
//
//  Why per bin is the right thing, measured not assumed:
//    * the kinematic acceptance inside a theta_CM bin is ~1.0 over Ex 0-16 (binacc.h computes it
//      exactly -- theta_CM and Ex fix theta_lab and KE completely), so NO acceptance correction
//      is needed inside a bin and acceptance cannot be blamed for the continuum shape;
//    * the continuum SHAPE differs between bins -- steeply falling at theta_CM 10-20, flat at
//      40-50 -- which is a cross-section angular dependence. The inclusive spectrum mixes them,
//      and that is why no single continuum shape ever fitted it.
//
//  Model per bin (the (p,p') recipe):
//    g.s. Gaussian + 0.740 5/2+ Gaussian, sigma FIXED to sigma_det
//    six barrier BWs, energies and Gamma FIXED, amplitudes free and bounded
//    continuum  floor_i + N_i * dens(Ex) * acc_i(Ex)
// =============================================================================================
namespace B2 {
   const double SN = 1.218, AMU = 931.494, HBARC = 197.327;
   const double R_FM = 1.4 * std::pow(15.0, 1.0 / 3.0), MU = AMU * 14.0 / 15.0;
   inline double P_ell(double Ed, int l)
   { if (Ed <= 0) return 0;
     double r = std::sqrt(2 * MU * Ed) / HBARC * R_FM;
     double F = r * ROOT::Math::sph_bessel(l, r), G = -r * ROOT::Math::sph_neumann(l, r);
     double d = F * F + G * G; return d > 0 ? r / d : 0; }
   inline double BW(double E, double Er, double Gr, int l)
   { double Ed = E - SN, Edr = Er - SN; if (Ed <= 0 || Edr <= 0) return 0;
     double Pr = P_ell(Edr, l); if (Pr <= 0) return 0;
     double h = 0.5 * Gr * P_ell(Ed, l) / Pr;
     return (h / TMath::Pi()) / ((E - Er) * (E - Er) + h * h); }
}
struct T2 {
   std::vector<double> y; double x0 = 0, dx = 1; int n = 0;
   void build(double Er, double Gr, int l, double sig, double lo, double hi, int N = 1200)
   { n = N; x0 = lo; dx = (hi - lo) / (N - 1); y.assign(N, 0.0);
     const double dE = 0.002, E0 = B2::SN + 1e-6, E1 = 23.0;
     int M = (int)((E1 - E0) / dE) + 1; std::vector<double> bw(M);
     for (int j = 0; j < M; j++) bw[j] = B2::BW(E0 + j * dE, Er, Gr, l);
     double cut = 5 * sig, nr = 1.0 / (std::sqrt(2 * TMath::Pi()) * sig);
     for (int i = 0; i < N; i++) { double x = x0 + i * dx, s = 0;
       int jl = std::max(0, (int)((x - cut - E0) / dE)), jh = std::min(M - 1, (int)((x + cut - E0) / dE));
       for (int j = jl; j <= jh; j++) { double u = (x - (E0 + j * dE)) / sig; s += bw[j] * std::exp(-0.5 * u * u); }
       y[i] = s * dE * nr; } }
   double Eval(double x) const
   { if (n < 2) return 0; double t = (x - x0) / dx; int i = (int)t;
     if (i < 0 || i >= n - 1) return 0; double f = t - i; return y[i] * (1 - f) + y[i + 1] * f; }
};
static const int NB = 6;
static const double EPP[NB] = {4.55, 5.05, 6.10, 6.80, 7.50, 8.40};
static const double GPP[NB] = {0.020, 0.020, 0.300, 0.300, 0.400, 0.500};
static const int    LPP[NB] = {1, 1, 2, 2, 2, 2};
static const double PPGS = 0.119, E52d = 0.740;
static T2 TT[NB];
static TGraph* GC = nullptr;      // dens(Ex)*acc_i(Ex) for the current bin
static double  FL = 0, SIGD = 0.4, GSMU = 0;
// response: bifurcated Gaussian CORE + exponential TAILS. Measured on the g.s. peak alone
// (lineshape.C): chi2/ndf 2.18 against 21.88 for a bifurcated Gaussian on ATTPCROOT, 7.00
// against 33.18 on Spyral. The plain Gaussian response was the single worst defect in every
// fit so far.
static double SGL = 0.4, SGR = 0.4, ALF_L = 1.5, ALF_R = 1.5;
// A LOW-ENERGY BACKGROUND, fitted ALONGSIDE the phase space rather than instead of it.
//   bkgLE(Ex) = B * erf((Ex - Sn)/(sqrt2 w)) * exp(-kappa * Ex)
// turning on at the neutron threshold and dying away, so it cannot touch the bound states and
// cannot masquerade as a resonance. Its origin is not established: part of it is a background
// localised at vz ~ 90-180 mm (an extra ~21 % of the events in that vz slice, in BOTH codes),
// but removing that by a vz cut leaves most of the excess, so there is more to it.
// p9 amplitude, p10 turn-on width, p11 decay constant.
static bool USE_LEBKG = false;
//  p0 gsAmp  p1 5/2+ amp  p2 continuum scale  p3..p8 BW amps  p9..p11 low-energy background
double fBin(double* x, double* p)
{
   double X = x[0];
   double v = FL + std::max(0.0, p[2]) * (GC ? std::max(0.0, GC->Eval(X)) : 0.0);
   if (USE_LEBKG) {
      double thr = 0.5 * (1 + std::erf((X - B2::SN) / (M_SQRT2 * std::max(0.05, p[10]))));
      v += std::max(0.0, p[9]) * thr * std::exp(-std::max(0.0, p[11]) * X);
   }
   v += p[0] * shapeGT(X, GSMU, SGL, SGR, ALF_L, ALF_R);
   v += p[1] * shapeGT(X, GSMU + E52d, SGL, SGR, ALF_L, ALF_R);
   for (int i = 0; i < NB; i++) v += std::max(0.0, p[3 + i]) * TT[i].Eval(X);
   return v;
}

static double BEST_NLL = -1;
struct BR { int bin; double cmLo, cmHi, N, nll, chi2ndf, cs, floor, bgB, bgW, bgK;
            double amp[NB], err[NB], gsA, e5A; TH1D* h; double A[NB]; };

// Reference amplitudes from the INCLUSIVE spectrum, used to seed and bound every bin.
// Without this the per-bin fit collapses: several bins drove all six amplitudes to zero and let
// the continuum absorb the resonances (Spyral bin 1 did exactly that).  This is the stabilisation
// the (p,p') analysis uses -- seed = reference x (bin counts / total), bound tightly above.
static double REF[NB] = {0}, REFTOT = 0;
BR fitBin(const Cfg& k, int ib, double exLo, double exHi, bool talk, bool useRef = true)
{
   Pars A = parsOf(k);
   Win w{14, 86, k.keLo, k.keHi};
   double cmLo = TCM0 + ib * TCMSTEP, cmHi = cmLo + TCMSTEP;
   TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
   TH1D* h = (TH1D*)f.Get(Form("%s_cm%d", k.nm, ib));
   BR o{}; o.bin = ib; o.cmLo = cmLo; o.cmHi = cmHi;
   if (!h) return o;
   h = (TH1D*)h->Clone(Form("c_%s_%d", k.nm, ib)); h->SetDirectory(nullptr);
   o.N = h->GetEntries(); o.h = h;
   if (o.N < 400) return o;

   SIGD = k.fw / 2.355;
   // response measured on THIS bin's own g.s. peak, then frozen for the rest of the fit
   {
      double pk0 = h->GetBinContent(h->GetMaximumBin());
      TF1 gt("gtb", f_GT, -2.0, 2.2, 8);
      gt.SetParameters(pk0, 0, SIGD, SIGD, 1.5, 1.5, 5, 0);
      gt.SetParLimits(1, -0.8, 0.8);
      gt.SetParLimits(2, 0.05, 2.0); gt.SetParLimits(3, 0.05, 2.0);
      gt.SetParLimits(4, 0.3, 6.0);  gt.SetParLimits(5, 0.3, 6.0);
      h->Fit(&gt, "RQ0NL");
      GSMU = gt.GetParameter(1); SGL = gt.GetParameter(2); SGR = gt.GetParameter(3);
      ALF_L = gt.GetParameter(4); ALF_R = gt.GetParameter(5);
      SIGD = SGR;
   }
   // floor: measured below Sn in THIS bin
   double fl = 0; int nf = 0;
   for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
      if (x > -2.5 && x < -1.5) { fl += h->GetBinContent(i); nf++; } }
   FL = nf ? fl / nf : 0; o.floor = FL;

   for (int i = 0; i < NB; i++) TT[i].build((EPP[i] - PPGS) + GSMU, GPP[i], LPP[i], SIGD, exLo - 1.0, 23.0);
   GC = new TGraph();
   for (int i = 0; i <= 460; i++) { double X = i * 0.05;
      GC->SetPoint(i, X, PS::dens(X, k.eb) * accBin(A, k.eb, X, cmLo, cmHi, w)); }

   BEST_NLL = -1;
   double pk = h->GetBinContent(h->GetMaximumBin());
   double c9 = GC->Eval(9.0);
   double lvl = 0; int nl = 0;
   for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
      if (x > 8 && x < 10) { lvl += h->GetBinContent(i); nl++; } }
   lvl = nl ? lvl / nl : 1;

   TF1* g = new TF1(Form("b_%s_%d_%d_%d", k.nm, ib, (int)(exLo * 10), (int)USE_LEBKG), fBin, exLo, exHi, 12);
   // If the window starts above the bound states, drop them: the g.s. lineshape is NOT Gaussian
   // (chi2/bin ~ 9 in its core) and letting it sit inside the window poisons the resonances.
   bool gsIn = true;   // the g.s. is ALWAYS in the model now -- its tail is what matters
   if (gsIn) { g->SetParameter(0, pk);  g->SetParLimits(0, 0, 1e9);
               g->SetParameter(1, pk * 0.08); g->SetParLimits(1, 0, 1e9); }
   else      { g->FixParameter(0, 0);   g->FixParameter(1, 0); }
   g->SetParameter(2, c9 > 0 ? std::max(0.0, lvl - FL) / c9 : 1.0); g->SetParLimits(2, 0, 1e7);
   double ampScale = (useRef && REFTOT > 0) ? o.N / REFTOT : 0;
   for (int i = 0; i < NB; i++) {
      double sc = TT[i].Eval((EPP[i] - PPGS) + GSMU);
      if (useRef && REFTOT > 0 && REF[i] > 0) {
         double seed = REF[i] * ampScale;
         g->SetParameter(3 + i, seed);
         g->SetParLimits(3 + i, 0, std::max(3.0 * seed, 0.4 * REF[i]));
      } else {
         g->SetParameter(3 + i, sc > 0 ? std::max(1.0, 0.15 * lvl) / sc : 1.0);
         g->SetParLimits(3 + i, 0, 1e9);
      }
   }
   if (!USE_LEBKG) { g->FixParameter(9, 0); g->FixParameter(10, 1); g->FixParameter(11, 0); }
   g->SetNpx(1500);
   TFitResultPtr r;
   if (!USE_LEBKG) { r = h->Fit(g, "RQ0NLS"); r = h->Fit(g, "RQ0NLS"); }
   else {
      // MULTI-START on the background: one seed put Spyral bin 1 in a minimum WORSE than the
      // model nested inside it (d(-2lnL) = -805, impossible). Try several, remember which SEED
      // won, then re-run from it -- restoring the parameters and re-fitting let the minimiser
      // wander back off, and fixing them instead zeroed every error.
      double best = 1e30; int bs = 0, is = 0;
      const double BS[3] = {0.05, 0.3, 1.0}, KS[3] = {0.10, 0.35, 0.9};
      auto seedBkg = [&](int idx) {
         // seed EVERY parameter, not just the background ones: p0/p1/p2 carried over from the
         // previous loop iteration, so re-running from the winning seed did not reproduce its
         // minimum and reported a likelihood that disagreed with the chi2.
         g->SetParameter(0, pk);        g->SetParLimits(0, 0, 1e9);
         g->SetParameter(1, pk * 0.08); g->SetParLimits(1, 0, 1e9);
         g->SetParameter(2, c9 > 0 ? std::max(0.0, lvl - FL) / c9 : 1.0); g->SetParLimits(2, 0, 1e7);
         g->SetParameter(9, BS[idx / 3] * lvl); g->SetParLimits(9, 0, 1e6);
         g->SetParameter(10, 1.0);              g->SetParLimits(10, 0.05, 3.0);
         g->SetParameter(11, KS[idx % 3]);      g->SetParLimits(11, 0.0, 2.0);
         for (int i = 0; i < NB; i++) if (useRef && REFTOT > 0 && REF[i] > 0)
            g->SetParameter(3 + i, REF[i] * o.N / REFTOT); };
      for (is = 0; is < 9; is++) {
         seedBkg(is);
         TFitResultPtr rr = h->Fit(g, "RQ0NLS");
         double v = 2 * rr->MinFcnValue();
         if (v < best) { best = v; bs = is; } }
      seedBkg(bs);
      r = h->Fit(g, "RQ0NLS");
      BEST_NLL = best;
   }

   double c2 = 0; int nd = 0;
   for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i), n = h->GetBinContent(i);
      if (x < exLo || x > exHi || n < 10) continue;
      double m = g->Eval(x); c2 += (n - m) * (n - m) / n; nd++; }
   int npar = (gsIn ? 9 : 7) + (USE_LEBKG ? 3 : 0);
   o.nll = (USE_LEBKG && BEST_NLL > 0) ? BEST_NLL : 2 * r->MinFcnValue();
   o.chi2ndf = c2 / std::max(1, nd - npar);
   o.cs = g->GetParameter(2); o.gsA = g->GetParameter(0); o.e5A = g->GetParameter(1);
   o.bgB = g->GetParameter(9); o.bgW = g->GetParameter(10); o.bgK = g->GetParameter(11);
   for (int i = 0; i < NB; i++) { o.A[i] = g->GetParameter(3 + i);
      double t = TT[i].Eval((EPP[i] - PPGS) + GSMU);
      o.amp[i] = o.A[i] * t; o.err[i] = g->GetParError(3 + i) * t; }
   if (talk) {
      printf("  bin %d  cm %2.0f-%2.0f  N %6.0f  chi2/ndf %6.2f  cont %.4g  floor %4.1f | ",
             ib, cmLo, cmHi, o.N, o.chi2ndf, o.cs, o.floor);
      for (int i = 0; i < NB; i++)
         printf(" BW%d %5.1f", i + 1, o.err[i] > 0 ? o.amp[i] / o.err[i] : 0.0);
      printf("\n");
   }
   return o;
}
// ADOPTED: Ex window [-2, 12] with the tailed response.  The scan over
// {-2..12, 2.5..12, 2.5..14, 3..12, 2.5..10, 2.5..16} makes the full window best once the g.s.
// tail is modelled (median chi2/ndf 2.77 / 2.02 against 10-18 for the truncated windows): cutting
// the g.s. out removes the constraint that pins the tail under the resonances.
void binfit()
{
   for (auto& k : CFG) {
      printf("\n######## %s : low-energy background ON vs OFF ########\n", k.nm);
      printf("  bin      N     chi2 OFF   chi2 ON    d(-2lnL)   B      w     kappa |  BW significances ON\n");
      for (int ib = 0; ib < 4; ib++) {
         USE_LEBKG = false;
         REFTOT = 0; for (int i = 0; i < NB; i++) REF[i] = 0;
         BR a = fitBin(k, ib, -2.0, 12.0, false, false);
         if (a.N < 400) continue;
         USE_LEBKG = true;
         BR b = fitBin(k, ib, -2.0, 12.0, false, false);
         printf("  %2d %7.0f %10.2f %9.2f %10.1f", ib, a.N, a.chi2ndf, b.chi2ndf, a.nll - b.nll);
         printf("  %6.1f %5.2f %6.3f |", b.bgB, b.bgW, b.bgK);
         for (int i = 0; i < NB; i++) printf(" %5.1f", b.err[i] > 0 ? b.amp[i] / b.err[i] : 0.0);
         printf("\n");
      }
      USE_LEBKG = false;
      // 1) inclusive reference: fit the sum of all theta_CM bins, free amplitudes
      REFTOT = 0; for (int i = 0; i < NB; i++) REF[i] = 0;
      {
         TFile fi("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
         TH1D* ha = (TH1D*)fi.Get(Form("%s_all", k.nm));
         double tot = ha ? ha->GetEntries() : 0;
         // reuse fitBin on a pseudo-bin covering the whole CM range by summing the bins
         BR acc{}; double sumA[NB] = {0};
         for (int ib = 0; ib < 6; ib++) {
            BR o = fitBin(k, ib, -2.0, 12.0, false, false);
            if (o.N < 400) continue;
            for (int i = 0; i < NB; i++) sumA[i] += o.A[i];
            acc.N += o.N;
         }
         for (int i = 0; i < NB; i++) REF[i] = sumA[i];
         REFTOT = acc.N;
         printf("\n[%s] reference amplitudes from the free per-bin pass, total N %.0f\n", k.nm, REFTOT);
      }
      printf("\n================ %s   per theta_CM bin, Ex [-2, 12] ================\n", k.nm);
      printf(" bin  cm      N    chi2  cont      |");
      for (int i = 0; i < NB; i++) printf("  BW%d(%.2f)", i + 1, EPP[i] - PPGS);
      printf("\n");
      double sumA[NB] = {0};
      for (int ib = 0; ib < 6; ib++) {
         BR o = fitBin(k, ib, -2.0, 12.0, false);
         if (o.N < 400) { printf(" %2d  %2.0f-%2.0f %6.0f   (too few)\n", ib, o.cmLo, o.cmHi, o.N); continue; }
         printf(" %2d  %2.0f-%2.0f %6.0f %6.2f %9.3g |", ib, o.cmLo, o.cmHi, o.N, o.chi2ndf, o.cs);
         for (int i = 0; i < NB; i++) {
            double sg = o.err[i] > 0 ? o.amp[i] / o.err[i] : 0.0;
            printf(" %6.1f%s", sg, sg > 5 ? "*" : " ");
            sumA[i] += o.amp[i];
         }
         printf("\n");
      }
      printf("  (* = above 5 sigma)   summed peak height:");
      for (int i = 0; i < NB; i++) printf(" %7.0f", sumA[i]);
      printf("\n");
   }
}
