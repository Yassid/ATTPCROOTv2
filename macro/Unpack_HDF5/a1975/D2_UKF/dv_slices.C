/// @file dv_slices.C
/// @brief Print E_x against theta_cm for several drift velocities side by side, so the flatness
///        can be judged by eye instead of through a zero crossing that keeps jumping.
///
/// The elastic 16C(d,d)16C channel has E_x = 0 at EVERY angle, so the column that is flattest
/// (and closest to zero) is the right dv. The statistic is the mean of E_x inside a +-1.5 MeV
/// window around zero, which is far steadier than the bin-quantised peak position the automated
/// scan in dt_dvscan.C used -- with only five or six slices in the trusted band, that peak
/// quantisation was enough to swing the fitted slope from -0.3 to +0.4 and send the reported
/// crossing anywhere between 1.00 and 1.42.
///
/// Below theta_cm ~ 30 the numbers are not usable: inverse kinematics puts those recoils at
/// theta_lab ~ 80 deg with very little energy, giving short near-transverse tracks. That is the
/// strong low-angle dependence seen in the viewer, and it is a reconstruction limit rather than
/// a scale error, so it must not drive the calibration.
///
///   root -b -q 'dv_slices.C(191)'

static double omS(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

void dv_slices(double Ebeam = 191.0, TString cache = "deuteron_kin_dd.root", double dvUsed = 1.15,
               double vzLo = 0, double vzHi = 500, double keMin = 3, double chi2max = 5)
{
   const double u = 931.49401;
   const double m1 = 16.0147013 * u, m2 = 2.0135532 * u, m3 = 2.0135532 * u, m4 = 16.0147013 * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   float ke, th, vz = 0, c2;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);

   std::vector<double> pT, pz, vzv;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (c2 > chi2max || ke <= 0) continue;
      double p = std::sqrt((ke + m3) * (ke + m3) - m3 * m3), a = th * TMath::DegToRad();
      pT.push_back(p * std::sin(a));
      pz.push_back(p * std::cos(a));
      vzv.push_back(vz);
   }

   const double dvs[] = {1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.40};
   const int nd = sizeof(dvs) / sizeof(double);
   const int aLo = 30, aHi = 60, aSt = 5, nsl = (aHi - aLo) / aSt;
   std::vector<std::vector<double>> val(nd, std::vector<double>(nsl, 0));
   std::vector<std::vector<long>> cnt(nd, std::vector<long>(nsl, 0));

   for (int j = 0; j < nd; ++j) {
      const double k = dvs[j] / dvUsed;
      for (size_t i = 0; i < pT.size(); ++i) {
         double z = k * vzv[i];
         if (z < vzLo || z > vzHi) continue;
         double qz = k * pz[i], qT = pT[i];
         double kep = std::sqrt(qT * qT + qz * qz + m3 * m3) - m3;
         if (kep < keMin) continue;
         double thp = std::atan2(qT, qz);
         double E1 = Ebeam + m1, E3 = kep + m3, E4 = E1 + m2 - E3;
         double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, uu = m2 * m2 + m3 * m3 - 2 * m2 * E3;
         double ar = (std::cos(thp) * omS(s, m1 * m1, m2 * m2) * omS(uu, m2 * m2, m3 * m3) -
                      (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) / (2 * m2 * m2) + s + uu - m2 * m2;
         if (ar < 0) continue;
         double m4x = std::sqrt(ar), ex = m4x - m4;
         if (std::fabs(ex) > 1.5) continue;
         double tt = m2 * m2 + m4x * m4x - 2 * m2 * E4;
         double g = (s * s + s * (2 * tt - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
                     (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) / (omS(s, m1 * m1, m2 * m2) * omS(s, m3 * m3, m4x * m4x));
         if (g < -1 || g > 1) continue;
         double tcm = (TMath::Pi() - std::acos(g)) * TMath::RadToDeg();
         int b = (int)((tcm - aLo) / aSt);
         if (b < 0 || b >= nsl) continue;
         val[j][b] += ex;
         cnt[j][b]++;
      }
   }

   printf("\n(d,d) elastic, Ebeam %.0f, vertex z %.0f-%.0f mm.  <E_x> in |E_x|<1.5, per theta_cm slice.\n", Ebeam, vzLo, vzHi);
   printf("E_x is 0 at every angle, so the flattest column nearest zero is the right dv.\n\n");
   printf("%-10s", "theta_cm");
   for (int j = 0; j < nd; ++j) printf("  %5.2f", dvs[j]);
   printf("\n");
   for (int b = 0; b < nsl; ++b) {
      printf("%-10s", Form("%d-%d", aLo + b * aSt, aLo + (b + 1) * aSt));
      for (int j = 0; j < nd; ++j) printf("  %5.2f", cnt[j][b] ? val[j][b] / cnt[j][b] : 0.);
      printf("\n");
   }
   printf("%-10s", "spread");
   for (int j = 0; j < nd; ++j) {
      double lo = 1e9, hi = -1e9;
      for (int b = 0; b < nsl; ++b)
         if (cnt[j][b] > 50) { double v = val[j][b] / cnt[j][b]; lo = std::min(lo, v); hi = std::max(hi, v); }
      printf("  %5.2f", hi > lo ? hi - lo : 0.);
   }
   printf("\n%-10s", "|mean|");
   for (int j = 0; j < nd; ++j) {
      double s = 0; long n = 0;
      for (int b = 0; b < nsl; ++b) if (cnt[j][b] > 50) { s += val[j][b]; n += cnt[j][b]; }
      printf("  %5.2f", n ? std::fabs(s / n) : 0.);
   }
   printf("\n");
}
