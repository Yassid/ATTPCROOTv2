#pragma once
#include "binfit.C"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// =============================================================================================
//  SIMULTANEOUS fit of all theta_CM bins.
//
//  WHY.  Fitting each bin independently leaves nothing to couple them, and the amplitudes come
//  out erratic rather than forming an angular distribution:
//      ATTPCROOT BW3 peak height per bin:  35.3,  0,  7.7,  0,  0
//      ATTPCROOT BW1:                      91.8, 18.1, 13.5, 5.0, 0      (bin-to-bin x5)
//  Alternating zeros with large values is a fit pathology, not physics. The worst offender is the
//  BACKGROUND: its turn-on width w and decay kappa were refitted in EVERY bin, 15 free parameters
//  describing one detector effect, free to trade against the resonances differently each time.
//
//  WHAT IS SHARED vs PER BIN
//    shared  : background turn-on width w, background decay kappa            (2 parameters)
//    per bin : background normalisation B, g.s. amplitude, 5/2+ amplitude,
//              continuum scale, and the six BW amplitudes                    (10 x nBins)
//    fixed   : resonance energies and widths; the response (sigma_L, sigma_R, alpha_L, alpha_R
//              and the g.s. centroid), measured on each bin's OWN g.s. peak beforehand.
//
//  The likelihood is the Poisson -2lnL summed over every Ex bin of every theta_CM bin.
// =============================================================================================
namespace SF {
   const int NBIN = 5;
   const int NPB  = 10;                    // parameters per bin
   struct BinData {
      TH1D *h = nullptr;
      double gsMu, sgL, sgR, alL, alR;     // response, fixed
      std::vector<T2> T;                   // BW templates, fixed
      TGraph *gc = nullptr;                // dens*acc continuum shape
      double floor = 0;
      bool in[6];
      double N = 0;
   };
   std::vector<BinData> B;
   double exLo = -2.0, exHi = 9.0;
   double turnOn = 2.5;

   // ---- SMOOTH ANGULAR DISTRIBUTIONS -------------------------------------------------------
   // Sharing the background was not enough: the amplitudes still came out as 35.3, 0, 8.5, 0, 0.
   // The zeros are the tell -- they sit ON the lower bound, because six overlapping resonances can
   // trade freely inside each bin and the minimiser picks whatever combination fits that bin.
   //
   // So parametrise each resonance's angular distribution instead of fitting 5 free amplitudes:
   //        A_i(theta) = exp(c_i0 + c_i1 * u + c_i2 * u^2),   u = cos(theta_CM)
   // Three parameters per resonance instead of five, smooth by construction, and STRICTLY
   // POSITIVE -- which removes the zero-boundary pathology altogether. An inelastic angular
   // distribution is a smooth function of theta_CM, so this is a physical constraint, not a
   // numerical trick.
   bool SMOOTH = false;
   std::vector<double> UCM;                 // cos(theta_CM) at each bin centre
   double ampOf(int b, int i, const double *p)
   {
      if (!SMOOTH) return std::max(0.0, p[2 + b * NPB + 4 + i]);
      const double *c = p + 2 + NBIN * 4 + i * 3;     // layout below
      double u = UCM[b];
      double e = c[0] + c[1] * u + c[2] * u * u;
      if (e > 50) e = 50;                              // keep exp() finite
      return std::exp(e);
   }
   double model(int b, double X, const double *p)
   {
      const BinData &d = B[b];
      const double *q = p + 2 + b * (SMOOTH ? 4 : NPB);
      double w = std::max(0.05, p[0]), kap = std::max(0.0, p[1]);
      double v = d.floor + std::max(0.0, q[2]) * (d.gc ? std::max(0.0, d.gc->Eval(X)) : 0.0);
      double thr = 0.5 * (1 + std::erf((X - turnOn) / (M_SQRT2 * w)));
      v += std::max(0.0, q[3]) * thr * std::exp(-kap * X);
      v += q[0] * shapeGT(X, d.gsMu, d.sgL, d.sgR, d.alL, d.alR);
      v += q[1] * shapeGT(X, d.gsMu + E52d, d.sgL, d.sgR, d.alL, d.alR);
      for (int i = 0; i < 6; i++)
         if (d.in[i]) v += ampOf(b, i, p) * d.T[i].Eval(X);
      return v;
   }
   double fcn(const double *p)
   {
      double s = 0;
      for (int b = 0; b < (int)B.size(); b++) {
         TH1D *h = B[b].h;
         for (int i = 1; i <= h->GetNbinsX(); i++) {
            double x = h->GetBinCenter(i);
            if (x < exLo || x > exHi) continue;
            double n = h->GetBinContent(i), m = model(b, x, p);
            if (m < 1e-9) m = 1e-9;
            s += 2.0 * (m - n + (n > 0 ? n * std::log(n / m) : 0.0));   // Poisson -2lnL
         }
      }
      return s;
   }
}

// Prepare one bin: measure the response on its own g.s. peak, build the BW templates and the
// continuum shape, and read the floor from below Sn. Everything here is FIXED for the fit.
void sfPrep(const Cfg &k, int ib)
{
   SF::BinData d;
   TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
   TH1D *h = (TH1D *)f.Get(Form("%s_cm%d", k.nm, ib));
   if (!h) { SF::B.push_back(d); return; }
   d.h = (TH1D *)h->Clone(Form("sf_%s_%d", k.nm, ib)); d.h->SetDirectory(nullptr);
   d.N = d.h->GetEntries();
   double s0 = k.fw / 2.355, pk = d.h->GetBinContent(d.h->GetMaximumBin());
   TF1 gt(Form("sfgt_%s_%d", k.nm, ib), f_GT, -2.0, 2.2, 8);
   gt.SetParameters(pk, 0, s0, s0, 1.5, 1.5, 5, 0);
   gt.SetParLimits(1, -0.8, 0.8); gt.SetParLimits(2, 0.05, 2.0); gt.SetParLimits(3, 0.05, 2.0);
   gt.SetParLimits(4, 0.3, 6.0);  gt.SetParLimits(5, 0.3, 6.0);
   d.h->Fit(&gt, "RQ0NL");
   d.gsMu = gt.GetParameter(1); d.sgL = gt.GetParameter(2); d.sgR = gt.GetParameter(3);
   d.alL = gt.GetParameter(4);  d.alR = gt.GetParameter(5);
   double fl = 0; int nf = 0;
   for (int i = 1; i <= d.h->GetNbinsX(); i++) { double x = d.h->GetBinCenter(i);
      if (x > -2.5 && x < -1.5) { fl += d.h->GetBinContent(i); nf++; } }
   d.floor = nf ? fl / nf : 0;
   d.T.resize(6);
   for (int i = 0; i < 6; i++) {
      double mu = (EPP[i] - PPGS) + d.gsMu;
      d.in[i] = (mu > SF::exLo + d.sgR) && (mu < SF::exHi - d.sgR);
      d.T[i].build(mu, GPP[i], LPP[i], d.sgR, SF::exLo - 1.0, 23.0);
   }
   Pars A = parsOf(k); Win w{14, 86, k.keLo, k.keHi};
   double cmLo = TCM0 + ib * TCMSTEP, cmHi = cmLo + TCMSTEP;
   d.gc = new TGraph();
   for (int i = 0; i <= 460; i++) { double X = i * 0.05;
      d.gc->SetPoint(i, X, PS::dens(X, k.eb) * accBin(A, k.eb, X, cmLo, cmHi, w)); }
   SF::B.push_back(d);
}

// Fit all bins at once.  smooth = false -> free per-bin amplitudes (5 per resonance);
// smooth = true -> A_i(theta) = exp(c0 + c1 u + c2 u^2) with u = cos(theta_CM) (3 per resonance).
// Parameter layout
//   [0] bkg_w   [1] bkg_kappa
//   free   : [2 + b*10 + 0..9]  gs, e52, cont, bkgB, 6 amplitudes      (per bin)
//   smooth : [2 + b*4  + 0..3]  gs, e52, cont, bkgB                    (per bin)
//            [2 + NBIN*4 + i*3 + 0..2]  c0, c1, c2                     (per resonance)
double sfRun(const Cfg &k, bool smooth, bool talk)
{
   SF::B.clear(); SF::UCM.clear();
   for (int ib = 0; ib < SF::NBIN; ib++) sfPrep(k, ib);
   for (int ib = 0; ib < SF::NBIN; ib++)
      SF::UCM.push_back(std::cos((TCM0 + (ib + 0.5) * TCMSTEP) * M_PI / 180.0));
   SF::SMOOTH = smooth;

   int nper = smooth ? 4 : SF::NPB;
   int npar = 2 + SF::NBIN * nper + (smooth ? 18 : 0);
   ROOT::Math::Minimizer *m = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
   m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000);
   m->SetTolerance(0.1); m->SetPrintLevel(0); m->SetStrategy(1);
   ROOT::Math::Functor fnc(&SF::fcn, npar);
   m->SetFunction(fnc);
   m->SetLimitedVariable(0, "bkg_w",     1.0, 0.02, 0.05, 3.0);
   m->SetLimitedVariable(1, "bkg_kappa", 0.3, 0.01, 0.0,  2.0);
   std::vector<double> lvlOf(SF::NBIN, 1);
   for (int b = 0; b < SF::NBIN; b++) {
      SF::BinData &d = SF::B[b]; if (!d.h) continue;
      double pk = d.h->GetBinContent(d.h->GetMaximumBin());
      double lvl = 0; int nl = 0;
      for (int i = 1; i <= d.h->GetNbinsX(); i++) { double x = d.h->GetBinCenter(i);
         if (x > 6 && x < 8) { lvl += d.h->GetBinContent(i); nl++; } }
      lvl = nl ? lvl / nl : 1; lvlOf[b] = lvl;
      int o = 2 + b * nper;
      m->SetLimitedVariable(o + 0, Form("gs_%d", b),   pk,        pk * 1e-3,  0, pk * 3);
      m->SetLimitedVariable(o + 1, Form("e52_%d", b),  pk * 0.05, pk * 1e-4,  0, pk);
      m->SetLimitedVariable(o + 2, Form("cont_%d", b), 1e-3,      1e-5,       0, 1e3);
      m->SetLimitedVariable(o + 3, Form("bkgB_%d", b), 0.5 * lvl, lvl * 1e-3, 0, lvl * 20);
      if (!smooth)
         for (int i = 0; i < 6; i++) {
            double sc = d.in[i] ? d.T[i].Eval((EPP[i] - PPGS) + d.gsMu) : 0;
            double seed = (sc > 0) ? 0.2 * lvl / sc : 0;
            if (!d.in[i]) { m->SetFixedVariable(o + 4 + i, Form("bw%d_%d", i + 1, b), 0); continue; }
            m->SetLimitedVariable(o + 4 + i, Form("bw%d_%d", i + 1, b), seed,
                                  std::max(1e-6, seed * 1e-2), 0, seed * 200 + 1e3);
         }
   }
   if (smooth)
      for (int i = 0; i < 6; i++) {
         int o = 2 + SF::NBIN * 4 + i * 3;
         // seed c0 from a typical amplitude, c1 and c2 at zero (flat in angle)
         double sc = SF::B[0].in[i] ? SF::B[0].T[i].Eval((EPP[i] - PPGS) + SF::B[0].gsMu) : 1;
         double a0 = (sc > 0) ? 0.2 * lvlOf[0] / sc : 1;
         m->SetLimitedVariable(o + 0, Form("c0_%d", i + 1), std::log(std::max(1e-6, a0)), 0.05, -30, 30);
         m->SetLimitedVariable(o + 1, Form("c1_%d", i + 1), 0.0, 0.05, -40, 40);
         m->SetLimitedVariable(o + 2, Form("c2_%d", i + 1), 0.0, 0.05, -40, 40);
      }
   m->Minimize();
   const double *p = m->X();
   double tot = 0; int ntot = 0;
   if (talk) {
      printf("\n================ %s   %s  (%d parameters) ================\n", k.nm,
             smooth ? "SMOOTH angular distributions" : "free per-bin amplitudes", npar);
      printf("   shared background:  w = %.3f   kappa = %.4f /MeV     -2lnL = %.1f\n",
             p[0], p[1], m->MinValue());
      printf("   bin  cm      N   chi2/ndf |");
      for (int i = 0; i < 6; i++) printf("  BW%d(%.2f)", i + 1, EPP[i] - PPGS);
      printf("\n");
      for (int b = 0; b < SF::NBIN; b++) {
         SF::BinData &d = SF::B[b]; if (!d.h) continue;
         double c2 = 0; int nd = 0;
         for (int i = 1; i <= d.h->GetNbinsX(); i++) { double x = d.h->GetBinCenter(i), n = d.h->GetBinContent(i);
            if (x < SF::exLo || x > SF::exHi || n < 10) continue;
            double mm = SF::model(b, x, p); c2 += (n - mm) * (n - mm) / n; nd++; }
         tot += c2; ntot += nd;
         printf("   %2d  %2.0f-%2.0f %6.0f  %6.2f  |", b, TCM0 + b * TCMSTEP, TCM0 + (b + 1) * TCMSTEP,
                d.N, c2 / std::max(1, nd - nper));
         for (int i = 0; i < 6; i++) {
            double ph = d.in[i] ? SF::ampOf(b, i, p) * d.T[i].Eval((EPP[i] - PPGS) + d.gsMu) : 0;
            printf(" %9.1f", ph); }
         printf("\n");
      }
      printf("   overall chi2/ndf = %.3f  (%d bins, %d par)\n", tot / std::max(1, ntot - npar), ntot, npar);
   }
   double v = m->MinValue(); delete m; return v;
}

void simfit()
{
   for (auto &k : CFG) {
      double a = sfRun(k, false, true);
      double b = sfRun(k, true,  true);
      printf("\n   >>> %s : free %.1f  vs  smooth %.1f   d(-2lnL) = %+.1f for 12 fewer parameters\n",
             k.nm, a, b, b - a);
   }
}
