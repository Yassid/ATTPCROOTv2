/// @file bragg_fit_dp.h
/// @brief Species-blind BRAGG-CURVE fit for tracks that stop in the a1975 D2 gas (prototype, 2026-11-01).
///        The charge profile along the track (centroid charges at their 3D path coordinate from the vertex end) is
///        fitted with the energy lost per bin by an ion of FREE mass m and initial kinetic energy KE0, integrated with
///        the velocity-scaled proton CATIMA table (Z = 1: stopping depends on velocity only) and converted to charge
///        with the per-run constant k (MeV per charge unit). Bins containing a SATURATED hit (raw amplitude above
///        the GET ceiling) are predicted but excluded from the chi2. TMinuit, plain FCN (cling constraints, see
///        eom_brho_dp.h). Needs eom_brho_dp.h included first.
namespace bragg {
struct Profile {
   std::vector<double> s, q;   // bin centre path coordinate [mm] from the vertex end, charge (pad-corrected)
   std::vector<double> a, b;   // bin edges [mm]
   std::vector<char> use;      // 1 fitted; 0 = saturated (measurement is a lower limit); 2 = ignored (e.g. beam spray near the axis)
   double k{0};                // MeV per charge unit
   const eom::Table *tp{nullptr};
};
/// merge consecutive bins until each covers at least minLen mm of path (4 mm centroids are far too noisy along a
/// 3-4 pad wide spiral: neighbouring bins swing from ~0 to 3x the trend, 2026-11-01); a merged bin is used only if
/// all its parts are usable. Also returns the Bragg-rise ratio: charge/mm of the last bin over the median bin.
inline Profile rebin(const Profile &P, double minLen, double &braggRatio)
{
   Profile R;
   R.k = P.k; R.tp = P.tp;
   const int n = P.s.size();
   int i = 0;
   while (i < n) {
      int j = i;
      double q = 0, qs = 0;
      bool use = true;
      bool ign = false;
      while (j < n && (P.b[j] - P.a[i] < minLen || j == i) && (P.use[j] == 2) == (P.use[i] == 2)) {
         q += P.q[j]; qs += P.q[j] * P.s[j]; use = use && P.use[j] == 1; ign = P.use[j] == 2; ++j;
      }
      // a short tail joins the previous bin
      if (j == n && !R.s.empty() && P.b[n - 1] - P.a[i] < 0.5 * minLen) {
         double wPrev = R.q.back();
         R.s.back() = (R.s.back() * wPrev + qs) / std::max(1e-9, wPrev + q);
         R.q.back() += q; R.b.back() = P.b[n - 1]; if (R.use.back() == 1 && !use) R.use.back() = 0;
         break;
      }
      R.a.push_back(P.a[i]); R.b.push_back(P.b[j - 1]); R.q.push_back(q);
      R.s.push_back(q > 0 ? qs / q : 0.5 * (P.a[i] + P.b[j - 1])); R.use.push_back(ign ? 2 : (use ? 1 : 0));
      i = j;
   }
   std::vector<double> dens;
   for (size_t k = 0; k < R.q.size(); ++k) if (R.use[k] != 2) dens.push_back(R.q[k] / std::max(1e-9, R.b[k] - R.a[k]));
   std::vector<double> sd(dens);
   std::sort(sd.begin(), sd.end());
   braggRatio = sd.size() < 2 ? 0 : dens.back() / std::max(1e-9, sd[sd.size() / 2]);
   return R;
}
/// energy deposited in each bin for (m, KE0); returns the stopping path coordinate
inline double predict(const Profile &P, double m, double ke0, std::vector<double> &dep, double ds = 0.2)
{
   const int n = P.s.size();
   dep.assign(n, 0);
   double ke = ke0, s = 0, sStop = -1;
   int j = 0;
   const double sEnd = P.b.back() + 50;
   while (s < sEnd && ke > 0) {
      double de = std::min(ke, eom::sVel(*P.tp, std::max(ke, 0.005), m) * ds);
      double sm = s + 0.5 * ds;
      while (j < n && sm > P.b[j]) ++j;
      if (j < n && sm >= P.a[j]) dep[j] += de;
      ke -= de;
      s += ds;
      if (ke <= 0) sStop = s;
   }
   return sStop;
}
struct Chi {
   const Profile *P;
   mutable std::vector<double> dep;
   double operator()(const double *p) const
   {
      const double m = p[0], ke0 = p[1];
      if (m <= 100 || ke0 <= 0.01) return 1e30;
      predict(*P, m, ke0, dep);
      double c2 = 0;
      for (size_t j = 0; j < P->s.size(); ++j) {
         double pred = dep[j] / P->k, meas = P->q[j];
         // a bin predicted BELOW a saturated measurement is fine (the measurement is a lower limit)
         if (P->use[j] == 2) continue;
         if (!P->use[j]) { if (pred < meas) c2 += (meas - pred) * (meas - pred) / (0.04 * meas * meas + 1e-9); continue; }
         double sig2 = 0.15 * 0.15 * std::max(pred, meas) * std::max(pred, meas) + 0.02 * 0.02 * meas * meas + 1.0;
         c2 += (meas - pred) * (meas - pred) / sig2;
      }
      return c2;
   }
};
Chi *gChi = nullptr;
inline void fcn(int &, double *, double &f, double *p, int) { f = (*gChi)(p); }
/// fit; start from (m0, ke0). Returns true on success; m, ke0 updated; chi2/ndf out.
inline bool fit(const Profile &P, double &m, double &ke0, double &chi2, int &ndf)
{
   Chi chi{&P};
   gChi = &chi;
   TMinuit mn(2);
   mn.SetPrintLevel(-1);
   mn.SetFCN(fcn);
   int ierr = 0;
   double arg[2] = {1, 0};
   mn.mnexcm("SET ERR", arg, 1, ierr);
   mn.mnparm(0, "m", m, 0.05 * m, 300, 8000, ierr);
   mn.mnparm(1, "ke0", ke0, 0.05 * ke0, 0.01, 200, ierr);
   arg[0] = 2000; arg[1] = 0.01;
   mn.mnexcm("MIGRAD", arg, 2, ierr);
   double e, edm, errdef;
   int nv, nx, ist;
   mn.GetParameter(0, m, e);
   mn.GetParameter(1, ke0, e);
   mn.mnstat(chi2, edm, errdef, nv, nx, ist);
   ndf = -2;
   for (char u : P.use) ndf += (u == 1);
   return ist >= 1;
}
} // namespace bragg
