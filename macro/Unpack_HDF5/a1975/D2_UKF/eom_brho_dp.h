/// @file eom_brho_dp.h
/// @brief Equation-of-motion Brho for the a1975 D2 PID (prototype, 2026-10-31). A track's radius shrinks as it
///        slows, so a single circle underestimates the VERTEX Brho of strongly decelerating ejectiles (16C(d,t)
///        sim tritons: -5..-7 % at Brho 0.45-0.60). Here the trajectory is integrated in the uniform field with
///        CATIMA energy loss in the D2 gas for a given mass hypothesis and fitted to the hits (TMinuit):
///          mode 0 (planar): xy only, theta fixed (from the AtHelix3D hit refit); params x0 y0 phi0 p0
///          mode 1 (3D):     x y z together; params x0 y0 z0 phi0 theta p0 (z residuals weighted by zW)
///        Curvature kappa = 0.299792458 B / p_perp [1/mm, p in MeV/c]; KE -= dE/dx * ds (3D path).
///        Residual = distance to the nearest curve sample, truncated at 8 mm (blob/noise hits).
///        Interpreted use only: load the libraries INSIDE the calling macro (cling crashes with lambdas /
///        ROOT::Math::Functor here when they are loaded with -e) -- hence TMinuit and a plain FCN.
#include "TMinuit.h"
namespace eom {
struct Table {
   double lo{0.01}, hi{400};
   std::vector<double> s;
   double dedx(double ke) const
   {
      if (ke <= lo) return s.front();
      double x = (std::log(ke) - std::log(lo)) / (std::log(hi) - std::log(lo)) * (s.size() - 1);
      int i = std::min<int>(x, s.size() - 2);
      double f = x - i;
      return s[i] * (1 - f) + s[i + 1] * f;
   }
};
/// species k = 0 p, 1 d, 2 t
inline double massMeV(int k) { return k == 0 ? 938.272 : k == 1 ? 1875.613 : 2808.921; }
inline void makeTables(Table tab[3], double density = 6.61e-5)
{
   const int A[3] = {1, 2, 3};
   const double mu[3] = {1.00727646, 2.01355321, 3.01550071};
   for (int k = 0; k < 3; ++k) {
      AtTools::AtELossCATIMA el(density);
      el.SetProjectile(A[k], 1, mu[k]);
      std::vector<std::tuple<int, int, int>> mat;
      mat.emplace_back(2, 1, 1);
      el.SetMaterial(mat);
      tab[k].s.clear();
      for (int i = 0; i < 500; ++i)
         tab[k].s.push_back(el.GetdEdx(std::exp(std::log(tab[k].lo) + (std::log(tab[k].hi) - std::log(tab[k].lo)) * i / 499.)));
   }
}
/// MASS-FREE stopping (2026-11-01): for Z = 1 ions the electronic stopping depends on the velocity only, so the
/// proton table at the same gamma serves any mass m: S(KE, m) = S_p((gamma - 1) m_p), gamma = 1 + KE/m.
inline double sVel(const Table &tp, double ke, double m)
{
   return tp.dedx((ke / m) * 938.272);
}
/// mean stopping power [MeV/mm] over the first Lseg mm of 3D path from the vertex for |p0| [MeV/c] and mass m
inline double meanS(const Table &tp, double p0, double m, double Lseg, double ds = 0.5)
{
   double ke = std::sqrt(p0 * p0 + m * m) - m, ke0 = ke, s = 0;
   while (s < Lseg && ke > 0.01) { ke -= sVel(tp, ke, m) * std::min(ds, Lseg - s); s += ds; }
   if (ke < 0) ke = 0;
   return (ke0 - ke) / Lseg;
}
/// the mass whose predicted mean stopping over Lseg equals Smeas (bisection; meanS rises with m at fixed p)
inline double solveMass(const Table &tp, double p0, double Lseg, double Smeas, double mLo = 300, double mHi = 8000)
{
   double fLo = meanS(tp, p0, mLo, Lseg) - Smeas, fHi = meanS(tp, p0, mHi, Lseg) - Smeas;
   if (fLo > 0) return mLo;
   if (fHi < 0) return mHi;
   for (int it = 0; it < 40; ++it) {
      double mm = std::sqrt(mLo * mHi), f = meanS(tp, p0, mm, Lseg) - Smeas;
      if (f > 0) mHi = mm; else mLo = mm;
   }
   return std::sqrt(mLo * mHi);
}
struct Chi {
   const std::vector<double> *hx, *hy, *hz;
   const Table *tab;
   double m, B, L, ds, zW, thFix;
   int sense, mode;
   bool velScaled{false}; // tab is the PROTON table, scaled by velocity to mass m
   mutable std::vector<double> cx, cy, cz;
   double operator()(const double *p) const
   {
      double x = p[0], y = p[1], z = 0, ph, th, pm;
      if (mode == 0) { ph = p[2]; pm = p[3]; th = thFix; }
      else { z = p[2]; ph = p[3]; th = p[4]; pm = p[5]; }
      const double sth = std::sin(th), cth = std::cos(th);
      if (sth < 0.02 || pm <= 0) return 1e30;
      double ke = std::sqrt(pm * pm + m * m) - m;
      // mode 0 steps ds in xy (path ds/sth), mode 1 steps ds in 3D
      const double dxy = mode == 0 ? ds : ds * sth, d3 = mode == 0 ? ds / sth : ds;
      int n = int(L / dxy) + 2;
      cx.resize(n); cy.resize(n); cz.resize(n);
      int k = 0;
      for (; k < n; ++k) {
         cx[k] = x; cy[k] = y; cz[k] = z;
         if (ke <= 0.01) { ++k; break; }
         double pp = std::sqrt(ke * ke + 2 * ke * m) * sth;
         ph += sense * 0.299792458 * B / pp * dxy;
         x += dxy * std::cos(ph); y += dxy * std::sin(ph); z += d3 * cth;
         ke -= (velScaled ? sVel(*tab, ke, m) : tab->dedx(ke)) * d3;
      }
      cx.resize(k); cy.resize(k); cz.resize(k);
      const int mm = k, nh = hx->size();
      double s2 = 0;
      int j = 0;
      for (int i = 0; i < nh; ++i) {
         int lo = std::max(0, j - 20), hi = std::min(mm - 1, j + 40), best = lo;
         double bd = 1e18;
         for (int q = lo; q <= hi; ++q) {
            double d = (cx[q] - (*hx)[i]) * (cx[q] - (*hx)[i]) + (cy[q] - (*hy)[i]) * (cy[q] - (*hy)[i]);
            if (mode == 1) d += zW * zW * (cz[q] - (*hz)[i]) * (cz[q] - (*hz)[i]);
            if (d < bd) { bd = d; best = q; }
         }
         j = best;
         s2 += std::min(bd, 64.0);
      }
      return s2;
   }
};
Chi *gChi = nullptr;
inline void fcn(int &, double *, double &f, double *p, int) { f = (*gChi)(p); }
/// hits ordered from the vertex end; start values from the end circle (xc,yc,R,sense) and theta.
/// Returns |p0| [MeV/c] (0 on failure); chi2 -> fmin; fitted theta -> thOut (mode 1, else theta).
inline double fitM(const std::vector<double> &hx, const std::vector<double> &hy, const std::vector<double> &hz,
                   double xc, double yc, double R, int sense, double theta, double mass, const Table &tp, int mode,
                   double &fmin, double &thOut, double B = 2.85, double zW = 1.0);
inline double fit(const std::vector<double> &hx, const std::vector<double> &hy, const std::vector<double> &hz,
                  double xc, double yc, double R, int sense, double theta, int k, const Table tab[3], int mode,
                  double &fmin, double &thOut, double B = 2.85, double zW = 1.0)
{
   return fitM(hx, hy, hz, xc, yc, R, sense, theta, -k - 1, tab[k], mode, fmin, thOut, B, zW);
}
/// mass >= 0: continuous mass on the velocity-scaled proton table tp; mass = -k-1: species table tp (k = 0 p, 1 d, 2 t)
inline double fitM(const std::vector<double> &hx, const std::vector<double> &hy, const std::vector<double> &hz,
                   double xc, double yc, double R, int sense, double theta, double mass, const Table &tp, int mode,
                   double &fmin, double &thOut, double B, double zW)
{
   const bool vs = mass >= 0;
   const double mUse = vs ? mass : massMeV(int(-mass - 1));
   const double sth = std::sin(theta);
   double a0 = std::atan2(hy[0] - yc, hx[0] - xc);
   double x0 = xc + R * std::cos(a0), y0 = yc + R * std::sin(a0), phi0 = a0 + sense * M_PI / 2;
   double p0 = 0.299792458 * B * R / sth;
   double ph = 0, prev = a0, phMax = 0;
   for (size_t i = 1; i < hx.size(); ++i) {
      double a = std::atan2(hy[i] - yc, hx[i] - xc), d = a - prev;
      while (d > M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      ph += d; prev = a; phMax = std::max(phMax, std::abs(ph));
   }
   const double Lxy = R * phMax + 20;
   Chi chi{&hx, &hy, &hz, &tp, mUse, B, mode == 0 ? Lxy : Lxy / sth, 1.0, zW, theta, sense, mode, vs};
   gChi = &chi;
   const int np = mode == 0 ? 4 : 6;
   TMinuit mn(np);
   mn.SetPrintLevel(-1);
   mn.SetFCN(fcn);
   int ierr = 0;
   double arg[2] = {1, 0};
   mn.mnexcm("SET ERR", arg, 1, ierr);
   int ip = 0;
   mn.mnparm(ip++, "x0", x0, 0.5, 0, 0, ierr);
   mn.mnparm(ip++, "y0", y0, 0.5, 0, 0, ierr);
   if (mode == 1) mn.mnparm(ip++, "z0", hz[0], 0.5, 0, 0, ierr);
   mn.mnparm(ip++, "phi0", phi0, 0.01, 0, 0, ierr);
   if (mode == 1) mn.mnparm(ip++, "theta", theta, 0.005, 0.02, M_PI - 0.02, ierr);
   mn.mnparm(ip++, "p0", p0, 0.02 * p0, 0, 0, ierr);
   arg[0] = 4000; arg[1] = 0.01;
   mn.mnexcm("MIGRAD", arg, 2, ierr);
   double edm, errdef;
   int nv, nx, ist;
   mn.mnstat(fmin, edm, errdef, nv, nx, ist);
   if (ist < 1) return 0;
   double v, e;
   mn.GetParameter(np - 1, v, e);
   thOut = theta;
   if (mode == 1) { double t2; mn.GetParameter(4, t2, e); thOut = t2; }
   return v;
}
} // namespace eom
