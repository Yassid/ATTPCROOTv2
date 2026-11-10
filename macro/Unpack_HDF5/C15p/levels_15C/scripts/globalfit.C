// ============================================================================================
//  GLOBAL fit of the 15C+d excitation-energy spectrum -- ATTPCROOT and Spyral, same raw data.
//
//  Components:  15C g.s. (1/2+)  +  0.740 MeV (5/2+)  +  state A  +  state B  +  background.
//
//  --------------------------------------------------------------------------------------------
//  1. Why TWO fit windows and not one simultaneous minimisation
//  --------------------------------------------------------------------------------------------
//  The g.s. carries ~9000 counts in its peak bin; state A carries ~300 in total. In a single
//  Poisson likelihood over the whole range the g.s. lineshape -- which is NOT Gaussian: sharp core,
//  longer right tail -- dominates the gradient and drags the weak states. Measured on this data:
//
//    one fit, free widths        A collapsed to the width floor and B blew up to FWHM 3.06
//                                (ATTPCROOT); both pinned at FWHM 0.283 (Spyral)
//    one fit, widths bounded     the background exponential slope railed for BOTH codes
//    erf x (exp + const) bkg     Spyral's A railed at 5.800 and B at 8.000 while the background
//                                scale hit its floor -- background and states swapped roles
//
//  The two regions are physically disjoint; the only thing linking them is the 5/2+ on the g.s.
//  right shoulder. So they are fitted in two windows SHARING ONE BACKGROUND, iterated to
//  self-consistency (the state window's background is frozen into the peak window, the peak
//  window's g.s. tail into the state window, 4 passes), then drawn as a single curve.
//
//  --------------------------------------------------------------------------------------------
//  2. The background is not known a priori, so it is a SYSTEMATIC, not a choice
//  --------------------------------------------------------------------------------------------
//  Three forms are fitted independently and the spread in the state energies is quoted as the
//  background systematic:
//
//    LIN    c0 + c1*Ex                                 empirical, monotonic, 2 parameters
//    PS     c + N*dens(Ex)*exp(-k*Ex)                  flat floor + three-body 15C+d -> d'+(14C+n)
//                                                      phase space (psmodel.h), threshold fixed
//                                                      at Sn, one acceptance factor, 3 parameters
//    THR    c + A*Phi((Ex-Sn)/sigma_res)*exp(-k*Ex)    flat floor + a continuum turning on at Sn
//                                                      smeared by that code's OWN measured g.s.
//                                                      sigma_R (not refitted), 3 parameters
//
//  None of the three can grow a Gaussian bump, so none can swap roles with a state. The result
//  (see the printout) is that state A and state B move by <= 0.06 MeV across all three forms in
//  both codes, while the 5/2+ amplitude ratio moves by 20-60 % -- so the states are robust and
//  the 5/2+ yield is NOT.
//
//  Pure threshold shapes cannot be the whole story: ATTPCROOT has ~70 counts/bin BELOW
//  Sn = 1.2180 MeV, which no 15C -> 14C+n shape produces, hence the flat floor in PS and THR.
//  Spyral has only ~7/bin there (it cuts vz < 400 mm and vertexes better), so the two codes
//  legitimately want different floors.
// ============================================================================================
#include "psmodel.h"

static const double E52 = 0.740;     // 15C 5/2+ excitation energy, MeV
static const double SNUC = 1.2180;   // 15C neutron separation energy, MeV -- FIXED everywhere
static const double SPLIT = 2.2;     // boundary between the peak and the state window, MeV

static int    MODE    = 0;           // 0 LIN, 1 PS, 2 THR
static double KB      = 191.25;      // beam KE feeding the phase-space shape, MeV
static double SIG_RES = 0.30;        // threshold smearing for THR = measured g.s. sigma_R, FIXED
static double P0, P1, P2;            // the shared background parameters

double bkgOf(double X, double a, double b, double c)
{
   if (MODE == 0) return std::max(0.0, a + b * X);
   if (MODE == 1) return std::max(0.0, a) + std::max(0.0, b) * PS::dens(X, KB) * std::exp(-c * X);
   double thr = 0.5 * (1.0 + std::erf((X - SNUC) / (M_SQRT2 * SIG_RES)));
   return std::max(0.0, a) + std::max(0.0, b) * thr * std::exp(-c * X);
}
double bkg(double X) { return bkgOf(X, P0, P1, P2); }

//  peak window:  p0 amp   p1 mu   p2 sigma_L   p3 sigma_R   p4 5/2+ / g.s. amplitude ratio
double fPeak(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   return bkg(X) + p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
        + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2));
}
static double GST[5];                // g.s. + 5/2+, frozen from the peak window
double gstail(double X)
{
   double s = (X < GST[1]) ? GST[2] : GST[3];
   return GST[0] * std::exp(-0.5 * std::pow((X - GST[1]) / s, 2))
        + GST[4] * GST[0] * std::exp(-0.5 * std::pow((X - GST[1] - E52) / GST[3], 2));
}
//  state window: p0,p1,p2 background   p3,p4,p5 state A (amp, mu, sigma)   p6,p7,p8 state B
double fStates(double* x, double* p)
{
   double X = x[0];
   return bkgOf(X, p[0], p[1], p[2]) + gstail(X)
        + p[3] * std::exp(-0.5 * std::pow((X - p[4]) / p[5], 2))
        + p[6] * std::exp(-0.5 * std::pow((X - p[7]) / p[8], 2));
}
//  assembled: p0-p4 peak, p5-p7 background, p8-p10 state A, p11-p13 state B
double fAll(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   return bkgOf(X, p[5], p[6], p[7])
        + p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
        + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2))
        + p[8] * std::exp(-0.5 * std::pow((X - p[9]) / p[10], 2))
        + p[11] * std::exp(-0.5 * std::pow((X - p[12]) / p[13], 2));
}

struct Out {
   double mu, dmu, fwL, fwR, r52, dr52, r52sys;
   double Amu, dAmu, Afw, dAfw, Bmu, dBmu, Bfw, dBfw;
   double floorFit, floorMeas, chi2ndf, c2core, c2shoulder, c2A, c2B, c2bkg;
   double par[14];
};

Out one(const char* nm, double kb, double gsFWHM, double lo, int mode, bool save)
{
   MODE = mode; KB = kb; SIG_RES = gsFWHM / 2.355;
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   double peak = h->GetMaximum(), s0 = gsFWHM / 2.355;

   auto avg = [&](double a, double b) { double s = 0; int n = 0;
      for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
         if (x > a && x < b) { s += h->GetBinContent(i); n++; } } return n ? s / n : 0.0; };
   double subThr = avg(-2.5, -1.5);               // measured, never given to the fit
   double yA = avg(8, 10), yB = avg(16, 20);
   if (mode == 0) { P1 = (yB - yA) / 9.0; P0 = yA - P1 * 9.0; P2 = 0; }
   else { P0 = std::max(0.0, subThr); P2 = (mode == 1) ? 0.15 : 0.04;
          double sh = (mode == 1) ? PS::dens(9.0, kb) : 1.0;
          P1 = std::max(1.0, (yA - P0) / (sh * std::exp(-P2 * 9.0))); }

   TF1* gp = new TF1(Form("P_%s_%d", nm, mode), fPeak, lo, SPLIT, 5);
   TF1* gs = new TF1(Form("S_%s_%d", nm, mode), fStates, SPLIT, 20.0, 9);
   gp->SetNpx(1000); gs->SetNpx(1000);

   for (int it = 0; it < 4; it++) {
      gp->SetParameter(0, peak); gp->SetParLimits(0, 0, 1e9);
      gp->SetParameter(1, 0);    gp->SetParLimits(1, -0.8, 0.8);
      gp->SetParameter(2, s0);   gp->SetParLimits(2, 0.05, 1.2);
      gp->SetParameter(3, s0);   gp->SetParLimits(3, 0.05, 1.2);
      gp->SetParameter(4, 0.08); gp->SetParLimits(4, 0, 0.8);
      h->Fit(gp, "RQ0NL");
      for (int i = 0; i < 5; i++) GST[i] = gp->GetParameter(i);

      gs->SetParameter(0, P0); gs->SetParLimits(0, 0, 1e4);
      if (mode == 0) { gs->SetParameter(1, P1); gs->SetParLimits(1, -50, 5);
                       gs->FixParameter(2, 0); }
      else           { gs->SetParameter(1, P1); gs->SetParLimits(1, 0, 1e7);
                       gs->SetParameter(2, P2); gs->SetParLimits(2, 0.0, 0.40); }
      gs->SetParameter(3, peak * 0.012); gs->SetParLimits(3, 0, 1e8);
      gs->SetParameter(4, 4.65);         gs->SetParLimits(4, 3.6, 5.8);
      gs->SetParameter(5, 0.60);         gs->SetParLimits(5, 0.10, 1.40);
      gs->SetParameter(6, peak * 0.012); gs->SetParLimits(6, 0, 1e8);
      gs->SetParameter(7, 6.88);         gs->SetParLimits(7, 6.0, 7.8);
      gs->SetParameter(8, 0.50);         gs->SetParLimits(8, 0.10, 1.40);
      h->Fit(gs, "RQ0NL");
      P0 = gs->GetParameter(0); P1 = gs->GetParameter(1); P2 = gs->GetParameter(2);
   }

   TF1* ga = new TF1(Form("A_%s_%d", nm, mode), fAll, lo, 20.0, 14);
   for (int i = 0; i < 5; i++) ga->SetParameter(i, gp->GetParameter(i));
   for (int i = 0; i < 9; i++) ga->SetParameter(5 + i, gs->GetParameter(i));
   ga->SetNpx(2000);

   auto chi = [&](double a, double b) { double c = 0; int n = 0;
      for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i), y = h->GetBinContent(i);
         if (x < a || x >= b || y < 10) continue; double m = ga->Eval(x); c += (y - m) * (y - m) / y; n++; }
      return std::make_pair(c, n); };
   auto tot = chi(lo, 20.0);

   Out o{};
   o.mu = gp->GetParameter(1); o.dmu = gp->GetParError(1);
   o.fwL = 2.355 * gp->GetParameter(2); o.fwR = 2.355 * gp->GetParameter(3);
   o.r52 = gp->GetParameter(4); o.dr52 = gp->GetParError(4); o.r52sys = 0;
   o.Amu = gs->GetParameter(4); o.dAmu = gs->GetParError(4);
   o.Afw = 2.355 * gs->GetParameter(5); o.dAfw = 2.355 * gs->GetParError(5);
   o.Bmu = gs->GetParameter(7); o.dBmu = gs->GetParError(7);
   o.Bfw = 2.355 * gs->GetParameter(8); o.dBfw = 2.355 * gs->GetParError(8);
   o.floorFit = (mode == 0) ? bkgOf(-1.9, P0, P1, P2) : P0;
   o.floorMeas = subThr;
   o.chi2ndf = tot.first / std::max(1, tot.second - 14);
   auto a1 = chi(lo, 1.5), a2 = chi(1.5, 3.5), a3 = chi(3.5, 5.8), a4 = chi(5.8, 8.0), a5 = chi(8.0, 20.0);
   o.c2core = a1.first / std::max(1, a1.second);  o.c2shoulder = a2.first / std::max(1, a2.second);
   o.c2A = a3.first / std::max(1, a3.second);     o.c2B = a4.first / std::max(1, a4.second);
   o.c2bkg = a5.first / std::max(1, a5.second);
   for (int i = 0; i < 14; i++) o.par[i] = ga->GetParameter(i);

   if (save) {
      TFile o2("/home/yassid/Data/C15d_statefit/globalfit.root", "UPDATE");
      h->Write(Form("h_%s", nm), TObject::kOverwrite);
      ga->Write(Form("all_%s", nm), TObject::kOverwrite);
      o2.Close();
   }
   return o;
}

void globalfit()
{
   { TFile o("/home/yassid/Data/C15d_statefit/globalfit.root", "RECREATE"); }
   const char* MN[3] = {"LIN", "PS ", "THR"};
   struct { const char* nm; double kb, fw, lo; } C[2] =
      {{"ATTPCROOT", 191.25, 1.012, -1.6}, {"Spyral", 195.00, 0.548, -1.2}};

   for (int c = 0; c < 2; c++) {
      printf("\n================ %s ================\n", C[c].nm);
      printf("bkg  chi2/ndf | g.s. mu        FWHM L/R   | 5/2+ ratio      | state A           FWHM"
             "        | state B           FWHM        | floor fit/meas\n");
      Out best{}; int bestm = 0; double bc = 1e30;
      Out all3[3];
      for (int m = 0; m < 3; m++) {
         Out o = one(C[c].nm, C[c].kb, C[c].fw, C[c].lo, m, false);
         all3[m] = o;
         printf("%s  %7.2f  | %+.4f(%.4f) %.3f %.3f | %.4f(%.4f) | %.3f(%.3f) %.3f(%.3f)"
                " | %.3f(%.3f) %.3f(%.3f) | %5.0f / %5.0f\n",
                MN[m], o.chi2ndf, o.mu, o.dmu, o.fwL, o.fwR, o.r52, o.dr52,
                o.Amu, o.dAmu, o.Afw, o.dAfw, o.Bmu, o.dBmu, o.Bfw, o.dBfw, o.floorFit, o.floorMeas);
         if (o.chi2ndf < bc) { bc = o.chi2ndf; best = o; bestm = m; }
      }
      double amin = 1e9, amax = -1e9, bmin = 1e9, bmax = -1e9, rmin = 1e9, rmax = -1e9;
      for (int m = 0; m < 3; m++) { amin = std::min(amin, all3[m].Amu); amax = std::max(amax, all3[m].Amu);
         bmin = std::min(bmin, all3[m].Bmu); bmax = std::max(bmax, all3[m].Bmu);
         rmin = std::min(rmin, all3[m].r52); rmax = std::max(rmax, all3[m].r52); }
      printf("  spread over the three background forms:  A %.3f MeV   B %.3f MeV   5/2+ ratio %.1f %%\n",
             amax - amin, bmax - bmin, 100.0 * (rmax - rmin) / (0.5 * (rmax + rmin)));
      printf("  ADOPTED: %s (lowest chi2/ndf %.2f)\n", MN[bestm], bc);
      printf("    g.s. 1/2+  %+.3f +- %.3f (stat) +- ---      FWHM %.3f / %.3f\n",
             best.mu, best.dmu, best.fwL, best.fwR);
      printf("    0.740 5/2+ amp/g.s. %.4f +- %.4f (stat)  +- %.4f (bkg)\n",
             best.r52, best.dr52, 0.5 * (rmax - rmin));
      printf("    state A    %.3f +- %.3f (stat) +- %.3f (bkg) MeV   FWHM %.3f +- %.3f\n",
             best.Amu, best.dAmu, 0.5 * (amax - amin), best.Afw, best.dAfw);
      printf("    state B    %.3f +- %.3f (stat) +- %.3f (bkg) MeV   FWHM %.3f +- %.3f\n",
             best.Bmu, best.dBmu, 0.5 * (bmax - bmin), best.Bfw, best.dBfw);
      printf("    chi2/bin by region:  g.s.core %.2f  shoulder(1.5-3.5) %.2f  A %.2f  B %.2f  bkg %.2f\n",
             best.c2core, best.c2shoulder, best.c2A, best.c2B, best.c2bkg);
      one(C[c].nm, C[c].kb, C[c].fw, C[c].lo, bestm, true);   // re-run the adopted one and store it
   }
}
