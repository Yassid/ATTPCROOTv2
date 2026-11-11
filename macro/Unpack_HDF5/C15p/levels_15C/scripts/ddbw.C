#include <Math/SpecFuncMathMore.h>
#include "psmodel.h"
// =============================================================================================
//  15C(d,d') fitted the way the 15C(p,p') analysis is fitted
//  (repo Yassid/a2091-15Cpp-analysis, macros/C15_pp_ana_AngBins_Barrier_6BW.C)
//
//  Model, mirroring theirs term for term:
//    - g.s. Gaussian (bound)                                        bifurcated, see note
//    - 0.740 MeV 5/2+ Gaussian (bound)
//    - SIX barrier-corrected Breit-Wigners, convolved with the detector Gaussian
//    - a 1n phase-space template
//  No empirical straight line anywhere.
//
//  Penetrability, identical to BarrierShape_demo.C / BarrierPhys:
//    rho   = k R,  k = sqrt(2 mu Edec)/hbar_c,  R = 1.4 * 15^(1/3) fm,  mu = 931.494*14/15
//    F_l   =  rho j_l(rho),   G_l = -rho y_l(rho),   P_l = rho/(F^2+G^2)
//    Gamma(E) = Gamma_r * P_l(E)/P_l(E_r)                     (Lane-Thomas)
//    BW(E)    = (Gamma(E)/2/pi) / [(E-E_r)^2 + (Gamma(E)/2)^2]
//
//  Resonance list and l-values are theirs (the 6-BW model):
//    BW1 4.55 l=1 | BW2 5.05 l=1 | BW3 6.10 l=2 | BW4 6.80 l=2 | BW5 7.50 l=2 | BW6 8.40 l=2
//  quoted in the (p,p') Ex scale, whose g.s. sits at +0.119. Each is therefore placed here at
//      E_here = (E_pp - 0.119) + (this code's own g.s. centroid)
//  i.e. every level is referred to its own code's ground state -- the only way the two codes
//  agree (see LEVELS_15C.md).
//
//  WHY the energies and widths are FIXED and only amplitudes float: the (d,d') resolution is
//  FWHM 0.91 (ATTPCROOT) and 0.52 (Spyral) against 0.51 for (p,p'), while BW1 and BW2 are 0.5 MeV
//  apart. Six free E_r and six free Gamma_r cannot be determined from this spectrum; fixing them
//  to the (p,p') values turns the fit into the question actually worth asking -- "is each (p,p')
//  resonance present in (d,d'), and how strong?" -- and makes every convolved BW a CONSTANT shape
//  that can be precomputed once, which is also what makes the fit affordable.
// =============================================================================================
namespace Bar {
   const double SN = 1.218, A_TRG = 14.0, AMU = 931.494, HBARC = 197.327;
   const double R_FM = 1.4 * std::pow(15.0, 1.0 / 3.0);
   const double MU = AMU * A_TRG / (A_TRG + 1.0);
   inline double rho_of(double Edec)
   { return Edec <= 0 ? 0.0 : std::sqrt(2.0 * MU * Edec) / HBARC * R_FM; }
   inline double P_ell(double Edec, int l)
   {
      double rho = rho_of(Edec);
      if (rho <= 0) return 0.0;
      double F = rho * ROOT::Math::sph_bessel(l, rho);
      double G = -rho * ROOT::Math::sph_neumann(l, rho);
      double d = F * F + G * G;
      return d > 0 ? rho / d : 0.0;
   }
   inline double BW(double E, double Er, double Gr, int l)
   {
      double Ed = E - SN, Edr = Er - SN;
      if (Ed <= 0 || Edr <= 0) return 0.0;
      double Pr = P_ell(Edr, l);
      if (Pr <= 0) return 0.0;
      double G_E = Gr * P_ell(Ed, l) / Pr, half = 0.5 * G_E;
      return (half / TMath::Pi()) / ((E - Er) * (E - Er) + half * half);
   }
}
// A convolved BW, precomputed on a grid (E_r, Gamma and sigma are fixed, so the shape is constant).
//
// ⚠ The convolution must be sampled on a grid fine compared to GAMMA, not to sigma. BW1 and BW2
// have Gamma_r = 20 keV while the detector sigma is 385 keV; sampling the Gaussian kernel at
// ~19 keV steps (201 points over +-5 sigma) under-samples a 20 keV Lorentzian and produced a
// visibly jagged template and an inflated chi2. Here the BW is cached once on a 2 keV grid and
// the Gaussian is integrated against THAT, which costs one Bessel pass per template and is both
// exact enough and fast.
struct Tmpl {
   std::vector<double> y; double x0 = 0, dx = 1; int n = 0;
   void build(double Er, double Gr, int l, double sig, double lo, double hi, int N = 1400)
   {
      n = N; x0 = lo; dx = (hi - lo) / (N - 1); y.assign(N, 0.0);
      const double dE = 0.002;                       // fine grid, 2 keV << Gamma_min = 20 keV
      const double E0 = Bar::SN + 1e-6, E1 = 23.0;
      int M = (int)((E1 - E0) / dE) + 1;
      std::vector<double> bw(M);
      for (int j = 0; j < M; j++) bw[j] = Bar::BW(E0 + j * dE, Er, Gr, l);
      double cut = 5.0 * sig, norm = 1.0 / (std::sqrt(2 * TMath::Pi()) * sig);
      for (int i = 0; i < N; i++) {
         double x = x0 + i * dx, sum = 0;
         int jlo = std::max(0, (int)((x - cut - E0) / dE));
         int jhi = std::min(M - 1, (int)((x + cut - E0) / dE));
         for (int j = jlo; j <= jhi; j++) {
            double u = (x - (E0 + j * dE)) / sig;
            sum += bw[j] * std::exp(-0.5 * u * u);
         }
         y[i] = sum * dE * norm;
      }
   }
   double Eval(double x) const
   {
      if (n < 2) return 0.0;
      double t = (x - x0) / dx; int i = (int)t;
      if (i < 0 || i >= n - 1) return 0.0;
      double f = t - i; return y[i] * (1 - f) + y[i + 1] * f;
   }
};

static const int NBW = 6;
static const double EPP[NBW] = {4.55, 5.05, 6.10, 6.80, 7.50, 8.40};   // (p,p') scale
static const double GPP[NBW] = {0.020, 0.020, 0.300, 0.300, 0.400, 0.500};
static const int    LPP[NBW] = {1, 1, 2, 2, 2, 2};
static const double PPGS = 0.119, E52 = 0.740;

static Tmpl TB[NBW];
static TH1D* HPS = nullptr;          // simulated + reconstructed 1n phase-space template
static const char* PSNAME = "A";     // "A" = d+14C+n, "B" = deuteron breakup (proton detected)
static bool  USEFLOOR = true;

//  p0 gsAmp p1 gsMu p2 sigL p3 sigR   p4 5/2+ ratio   p5 PS scale   p6 floor   p7..p12 BW amps
double fModel(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   double v = p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
            + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2));
   v += std::max(0.0, p[5]) * (HPS ? HPS->Interpolate(X) : 0.0);
   if (USEFLOOR) v += std::max(0.0, p[6]);
   for (int i = 0; i < NBW; i++) v += std::max(0.0, p[7 + i]) * TB[i].Eval(X);
   return v;
}

struct Res { double nll, chi2ndf; double amp[NBW], err[NBW], mu[NBW]; double gsMu, fwL, fwR, r52, ps, fl; };

Res fit(const char* nm, double fw, double lo, double shift, bool floorOn, bool talk)
{
   USEFLOOR = floorOn;
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   TFile fp("/home/yassid/Data/C15d_statefit/ps_alt.root");
   HPS = (TH1D*)fp.Get(PSNAME); HPS->SetDirectory(nullptr);

   double peak = h->GetBinContent(h->GetMaximumBin()), s0 = fw / 2.355;
   auto avg = [&](double a, double b) { double s = 0; int n = 0;
      for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
         if (x > a && x < b) { s += h->GetBinContent(i); n++; } } return n ? s / n : 0.0; };
   double subThr = avg(-2.5, -1.5), yA = avg(8, 10);

   // find the g.s. centroid first -- every resonance is placed relative to it
   TF1 gq("gq", "gaus(0)+pol1(3)", lo, 1.6);
   gq.SetParameters(peak, 0, s0, 10, 0); gq.SetParLimits(2, 0.05, 1.5);
   h->Fit(&gq, "RQ0N");
   double gsMu = gq.GetParameter(1);

   // `shift` is used as a SCALE on the excitation energy, not an offset: an offset just lets
   // BW5 slide into BW4's place and reshuffles the amplitudes, whereas a scale is what an
   // actual beam-energy / kinematic-scale error does.  E_here = gs + scale*(E_pp - pp gs).
   double SC = 1.0 + shift;
   for (int i = 0; i < NBW; i++)
      TB[i].build(SC * (EPP[i] - PPGS) + gsMu, GPP[i], LPP[i], s0, lo - 1.0, 23.0);

   TF1* g = new TF1(Form("bw_%s_%d", nm, (int)(shift * 1000)), fModel, lo, 20.0, 13);
   g->SetParameter(0, peak);    g->SetParLimits(0, 0, 1e9);
   g->SetParameter(1, gsMu);    g->SetParLimits(1, gsMu - 0.4, gsMu + 0.4);
   g->SetParameter(2, s0);      g->SetParLimits(2, 0.05, 1.2);
   g->SetParameter(3, s0);      g->SetParLimits(3, 0.05, 1.2);
   g->SetParameter(4, 0.08);    g->SetParLimits(4, 0, 0.8);
   g->SetParameter(5, std::max(1.0, yA / std::max(1e-9, HPS->Interpolate(9.0))));
   g->SetParLimits(5, 0, 1e6);
   if (floorOn) { g->SetParameter(6, subThr); g->SetParLimits(6, 0, 1e4); }
   else g->FixParameter(6, 0);
   for (int i = 0; i < NBW; i++) {
      double sc = TB[i].Eval(SC * (EPP[i] - PPGS) + gsMu);
      g->SetParameter(7 + i, sc > 0 ? 40.0 / sc : 1.0);
      g->SetParLimits(7 + i, 0, 1e9);
   }
   g->SetNpx(2000);
   TFitResultPtr r = h->Fit(g, "RQ0NLS"); r = h->Fit(g, "RQ0NLS");

   double c2 = 0; int nd = 0;
   for (int i = 1; i <= h->GetNbinsX(); i++) {
      double x = h->GetBinCenter(i), n = h->GetBinContent(i);
      if (x < lo || x > 20.0 || n < 10) continue;
      double m = g->Eval(x); c2 += (n - m) * (n - m) / n; nd++;
   }
   Res o{}; o.nll = 2 * r->MinFcnValue(); o.chi2ndf = c2 / std::max(1, nd - 13);
   o.gsMu = g->GetParameter(1); o.fwL = 2.355 * g->GetParameter(2);
   o.fwR = 2.355 * g->GetParameter(3); o.r52 = g->GetParameter(4);
   o.ps = g->GetParameter(5); o.fl = floorOn ? g->GetParameter(6) : 0;
   for (int i = 0; i < NBW; i++) {
      // peak height each BW contributes, which is what "significance" should be measured on
      double mu = SC * (EPP[i] - PPGS) + gsMu;
      o.mu[i] = mu;
      o.amp[i] = g->GetParameter(7 + i) * TB[i].Eval(mu);
      o.err[i] = g->GetParError(7 + i) * TB[i].Eval(mu);
   }
   if (talk) {
      printf("\n== %-10s  Ex scale %.4f  floor %s   -2lnL %.1f   chi2/ndf %.2f\n",
             nm, SC, floorOn ? "on" : "OFF", o.nll, o.chi2ndf);
      printf("   g.s. mu %+.3f  FWHM %.3f/%.3f   5/2+ ratio %.4f   PS scale %.3g   floor %.1f (measured %.1f)\n",
             o.gsMu, o.fwL, o.fwR, o.r52, o.ps, o.fl, subThr);
      printf("   BW   l   E_pp    E here    peak height +- err    significance\n");
      for (int i = 0; i < NBW; i++)
         printf("   BW%d  %d  %5.2f   %6.3f   %9.1f +- %5.1f    %6.1f sigma%s\n",
                i + 1, LPP[i], EPP[i], o.mu[i], o.amp[i], o.err[i],
                o.err[i] > 0 ? o.amp[i] / o.err[i] : 0.0,
                (o.err[i] > 0 && o.amp[i] / o.err[i] > 5) ? "  <<" : "");
   }
   return o;
}

void ddbw()
{
   struct { const char* nm; double fw, lo; } C[2] = {{"ATTPCROOT", 1.012, -1.6}, {"Spyral", 0.548, -1.2}};
   // Does the continuum want to be PHASE SPACE at all?  Four combinations, same resonances.
   printf("\n######## continuum: phase-space template only, vs template + flat floor ########\n");
   printf("%-10s  %-28s  %9s  %8s  %10s  %8s\n", "code", "continuum", "-2lnL", "chi2/ndf", "PS scale", "floor");
   for (auto& c : C) {
      for (const char* t : {"A", "B"}) {
         for (int fl = 0; fl < 2; fl++) {
            PSNAME = t;
            Res r = fit(c.nm, c.fw, c.lo, 0.0, fl == 1, false);
            printf("%-10s  %-28s  %9.1f  %8.2f  %10.4g  %8.1f\n", c.nm,
                   Form("%s (%s)%s", t, strcmp(t,"A")==0 ? "d+14C+n" : "breakup", fl ? " + floor" : " ONLY"),
                   r.nll, r.chi2ndf, r.ps, r.fl);
         }
      }
   }
   PSNAME = "A";
   printf("\n######## does (d,d') want the (p,p') energy SCALE? ########\n");
   for (auto& c : C) {
      double bs = 0, bl = 1e30;
      printf("%-10s  scale:", c.nm);
      for (double sc = -0.10; sc <= 0.1001; sc += 0.02) {
         Res r = fit(c.nm, c.fw, c.lo, sc, true, false);
         printf("  %.2f:%.0f", 1 + sc, r.nll);
         if (r.nll < bl) { bl = r.nll; bs = sc; }
      }
      printf("   best %.3f\n", 1 + bs);
      fit(c.nm, c.fw, c.lo, bs, true, true);
   }
}
