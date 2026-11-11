#include <Math/SpecFuncMathMore.h>
#include "psmodel.h"
// =============================================================================================
//  15C(d,d') -- the fit with the CORRECTED continuum.
//
//  Everything before this used a continuum template made with TGenPhaseSpace as an uncorrelated
//  three-body d + 14C + n final state. That is the wrong process: it throws the deuteron FORWARD
//  (0.1 % of its weight at theta 70-80, where 35 % of the data sits), so every fit discarded it.
//
//  Here the continuum is the SEQUENTIAL process actually at work,
//      15C + d -> d' + 15C*(M),   then   15C* -> 14C + n,
//  for which the deuteron has TWO-BODY kinematics against a recoil of mass M = M(15C) + Ex:
//
//      continuum(Ex) = floor + N * dens(Ex) * acc(Ex)
//
//        dens(Ex) = p3*(M) q(M)      the two-body/invariant-mass phase space of psmodel.h
//        acc(Ex)                     fraction of CM solid angle whose deuteron lands inside this
//                                    code's own (theta, KE) window -- pure kinematics, computed
//                                    here, no free parameter
//        floor                       FIXED to that code's MEASURED counts/bin below S_n, where
//                                    no 15C decay can contribute. Not a fit parameter.
//
//  Resonances: the six of the 15C(p,p') analysis, barrier-corrected BW (Lane-Thomas
//  penetrability) convolved with the detector resolution, energies and widths fixed.
// =============================================================================================
namespace Bar {
   const double SN = 1.218, A_TRG = 14.0, AMU = 931.494, HBARC = 197.327;
   const double R_FM = 1.4 * std::pow(15.0, 1.0 / 3.0);
   const double MU = AMU * A_TRG / (A_TRG + 1.0);
   inline double P_ell(double Edec, int l)
   {
      if (Edec <= 0) return 0.0;
      double rho = std::sqrt(2.0 * MU * Edec) / HBARC * R_FM;
      double F = rho * ROOT::Math::sph_bessel(l, rho), G = -rho * ROOT::Math::sph_neumann(l, rho);
      double d = F * F + G * G; return d > 0 ? rho / d : 0.0;
   }
   inline double BW(double E, double Er, double Gr, int l)
   {
      double Ed = E - SN, Edr = Er - SN;
      if (Ed <= 0 || Edr <= 0) return 0.0;
      double Pr = P_ell(Edr, l); if (Pr <= 0) return 0.0;
      double half = 0.5 * Gr * P_ell(Ed, l) / Pr;
      return (half / TMath::Pi()) / ((E - Er) * (E - Er) + half * half);
   }
}
struct Tmpl {
   std::vector<double> y; double x0 = 0, dx = 1; int n = 0;
   void build(double Er, double Gr, int l, double sig, double lo, double hi, int N = 1400)
   {
      n = N; x0 = lo; dx = (hi - lo) / (N - 1); y.assign(N, 0.0);
      const double dE = 0.002, E0 = Bar::SN + 1e-6, E1 = 23.0;
      int M = (int)((E1 - E0) / dE) + 1; std::vector<double> bw(M);
      for (int j = 0; j < M; j++) bw[j] = Bar::BW(E0 + j * dE, Er, Gr, l);
      double cut = 5 * sig, nrm = 1.0 / (std::sqrt(2 * TMath::Pi()) * sig);
      for (int i = 0; i < N; i++) {
         double x = x0 + i * dx, s = 0;
         int jl = std::max(0, (int)((x - cut - E0) / dE)), jh = std::min(M - 1, (int)((x + cut - E0) / dE));
         for (int j = jl; j <= jh; j++) { double u = (x - (E0 + j * dE)) / sig; s += bw[j] * std::exp(-0.5 * u * u); }
         y[i] = s * dE * nrm;
      }
   }
   double Eval(double x) const
   { if (n < 2) return 0; double t = (x - x0) / dx; int i = (int)t;
     if (i < 0 || i >= n - 1) return 0; double f = t - i; return y[i] * (1 - f) + y[i + 1] * f; }
};
// --- kinematic acceptance of the sequential channel, per code ---------------------------------
static double TH_LO = 14, TH_HI = 86, KE_LO = 0.5, KE_HI = 18.0;
double acc_of(double Ex, double Kb)
{
   double M4 = PS::M15C + Ex, M3 = PS::MD;
   double Eb = Kb + PS::M15C, pb = std::sqrt(Kb * Kb + 2 * Kb * PS::M15C);
   TLorentzVector W(0, 0, pb, Eb + PS::MD);
   double s = W.M2(), rs = std::sqrt(s);
   if (rs <= M3 + M4) return 0.0;
   double pcm = PS::omega(s, M3 * M3, M4 * M4) / (2 * rs), E3 = std::sqrt(pcm * pcm + M3 * M3);
   TVector3 b = W.BoostVector();
   const int N = 1500; double a = 0;
   for (int i = 0; i < N; i++) {
      double c = -1 + 2.0 * (i + 0.5) / N, st = std::sqrt(std::max(0.0, 1 - c * c));
      TLorentzVector d(pcm * st, 0, pcm * c, E3); d.Boost(b);
      double th = d.Theta() * 180 / M_PI, ke = d.E() - M3;
      if (th >= TH_LO && th <= TH_HI && ke >= KE_LO && ke <= KE_HI) a += 1;
   }
   return a / N;
}
static const int NBW = 6;
static const double EPP[NBW] = {4.55, 5.05, 6.10, 6.80, 7.50, 8.40};
static const double GPP[NBW] = {0.020, 0.020, 0.300, 0.300, 0.400, 0.500};
static const int    LPP[NBW] = {1, 1, 2, 2, 2, 2};
static const double PPGS = 0.119, E52 = 0.740;
static Tmpl TB[NBW];
static TGraph* GACC = nullptr;        // dens(Ex)*acc(Ex), precomputed
static double  FLOOR = 0;             // FIXED, measured below Sn
//  p0 gsAmp p1 gsMu p2 sigL p3 sigR  p4 5/2+ ratio  p5 continuum scale  p6..p11 BW amps
double fMod(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   double v = FLOOR + std::max(0.0, p[5]) * (GACC ? std::max(0.0, GACC->Eval(X)) : 0.0);
   v += p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2));
   v += p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2));
   for (int i = 0; i < NBW; i++) v += std::max(0.0, p[6 + i]) * TB[i].Eval(X);
   return v;
}

struct R { double nll, chi2ndf, gsMu, fwL, fwR, r52, cs, floor, sub;
           double amp[NBW], err[NBW], mu[NBW]; TH1D* h; double gsAmp, A[NBW]; };

R fitFinal(const char* nm, double kb, double fw, double lo,
           double thLo, double thHi, double keLo, double keHi, bool talk)
{
   TH_LO = thLo; TH_HI = thHi; KE_LO = keLo; KE_HI = keHi;
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   double peak = h->GetBinContent(h->GetMaximumBin()), s0 = fw / 2.355;
   auto avg = [&](double a, double b) { double s = 0; int n = 0;
      for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
         if (x > a && x < b) { s += h->GetBinContent(i); n++; } } return n ? s / n : 0.0; };
   double sub = avg(-2.5, -1.5);
   FLOOR = sub;                                   // FIXED, not fitted

   // continuum shape: dens * acc, precomputed
   GACC = new TGraph();
   for (int i = 0; i <= 460; i++) { double X = i * 0.05; GACC->SetPoint(i, X, PS::dens(X, kb) * acc_of(X, kb)); }

   TF1 gq("gq", "gaus(0)+pol1(3)", lo, 1.6);
   gq.SetParameters(peak, 0, s0, 10, 0); gq.SetParLimits(2, 0.05, 1.5);
   h->Fit(&gq, "RQ0N");
   double gsMu = gq.GetParameter(1);
   for (int i = 0; i < NBW; i++) TB[i].build((EPP[i] - PPGS) + gsMu, GPP[i], LPP[i], s0, lo - 1.0, 23.0);

   TF1* g = new TF1(Form("fin_%s", nm), fMod, lo, 20.0, 12);
   g->SetParameter(0, peak);  g->SetParLimits(0, 0, 1e9);
   g->SetParameter(1, gsMu);  g->SetParLimits(1, gsMu - 0.4, gsMu + 0.4);
   g->SetParameter(2, s0);    g->SetParLimits(2, 0.05, 1.2);
   g->SetParameter(3, s0);    g->SetParLimits(3, 0.05, 1.2);
   g->SetParameter(4, 0.08);  g->SetParLimits(4, 0, 0.8);
   double c9 = GACC->Eval(9.0);
   g->SetParameter(5, c9 > 0 ? std::max(0.0, avg(8, 10) - FLOOR) / c9 : 1.0); g->SetParLimits(5, 0, 1e6);
   for (int i = 0; i < NBW; i++) { double sc = TB[i].Eval((EPP[i] - PPGS) + gsMu);
      g->SetParameter(6 + i, sc > 0 ? 40.0 / sc : 1.0); g->SetParLimits(6 + i, 0, 1e9); }
   g->SetNpx(2000);
   TFitResultPtr r = h->Fit(g, "RQ0NLS"); r = h->Fit(g, "RQ0NLS");

   double c2 = 0; int nd = 0;
   for (int i = 1; i <= h->GetNbinsX(); i++) {
      double x = h->GetBinCenter(i), n = h->GetBinContent(i);
      if (x < lo || x > 20.0 || n < 10) continue;
      double m = g->Eval(x); c2 += (n - m) * (n - m) / n; nd++;
   }
   R o{}; o.nll = 2 * r->MinFcnValue(); o.chi2ndf = c2 / std::max(1, nd - 12);
   o.gsMu = g->GetParameter(1); o.fwL = 2.355 * g->GetParameter(2); o.fwR = 2.355 * g->GetParameter(3);
   o.r52 = g->GetParameter(4); o.cs = g->GetParameter(5); o.floor = FLOOR; o.sub = sub;
   o.h = h; o.gsAmp = g->GetParameter(0);
   for (int i = 0; i < NBW; i++) {
      double mu = (EPP[i] - PPGS) + gsMu; o.mu[i] = mu; o.A[i] = g->GetParameter(6 + i);
      o.amp[i] = g->GetParameter(6 + i) * TB[i].Eval(mu);
      o.err[i] = g->GetParError(6 + i) * TB[i].Eval(mu);
   }
   if (talk) {
      printf("\n== %-10s  continuum = floor(FIXED %.1f) + %.3g x dens(Ex)*acc(Ex)   -2lnL %.1f  chi2/ndf %.2f\n",
             nm, o.floor, o.cs, o.nll, o.chi2ndf);
      printf("   g.s. mu %+.3f  FWHM %.3f/%.3f   5/2+ ratio %.4f\n", o.gsMu, o.fwL, o.fwR, o.r52);
      printf("   continuum level: %.0f at Ex 2, %.0f at 9, %.0f at 18   (data %.0f / %.0f / %.0f)\n",
             FLOOR + o.cs * GACC->Eval(2), FLOOR + o.cs * GACC->Eval(9), FLOOR + o.cs * GACC->Eval(18),
             avg(1.5, 2.5), avg(8.5, 9.5), avg(17.5, 18.5));
      printf("   BW   l   E_x      peak height +- err    significance\n");
      for (int i = 0; i < NBW; i++)
         printf("   BW%d  %d  %5.2f   %9.1f +- %5.1f    %6.1f sigma%s\n", i + 1, LPP[i], o.mu[i],
                o.amp[i], o.err[i], o.err[i] > 0 ? o.amp[i] / o.err[i] : 0.0,
                (o.err[i] > 0 && o.amp[i] / o.err[i] > 5) ? "  <<" : "");
   }
   return o;
}
void ddfinal()
{
   fitFinal("ATTPCROOT", 191.25, 1.012, -1.6, 14, 86, 0.5, 18.0, true);
   fitFinal("Spyral",    195.00, 0.548, -1.2, 14, 86, 0.0, 20.0, true);
}
