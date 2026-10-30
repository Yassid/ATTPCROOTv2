/// @file dir3d_score.C
/// @brief Score every forward/backward vote of dir3d_bench.C against a referee that knows the vertex end.
///
/// Referee, per source:
///   sim (src 1): MC truth -- the track end nearest the generated vertex (3D, both z conventions
///                tried and the better one reported, not assumed).
///   H2  (src 0, mode 1): the 16C(p,p) elastic locus -- the end hypothesis that puts the track on
///                |Ex| < 2 when the other does not.
/// Votes (each names the VERTEX END, 0 or 1, or abstains):
///   rho    end nearer the beam axis
///   dca    end whose local helix passes closer to the axis (abstain if |diff| < dcaMin)
///   sense  end whose rotation sense, travelling inward, matches the CALIBRATED expectation
///   radius end with the larger local radius (decelerating spiral)
///   bragg  end with the lower edge charge density
///   spyral AtSpyralPID direction mapped to an end (forward = vertex at lowest z_digi)
///   atpra  GeoTheta > 90 -> vertex at lowest z_digi (AtGenfitter backwardSeed rule)
///   combo  dca if decisive, else sense
///
///   root -l -b -q 'pid/dir3d_score.C("/mnt/f/a1975/dir3d/d3_sim_dp.root",1)'
///   root -l -b -q 'pid/dir3d_score.C("/mnt/f/a1975/dir3d/d3_h2_0106.root",0)'

namespace {
double exPP(double Eb, double th, double ke)
{
   const double kU = 931.49410242, kMe = 0.51099895, kMp = 938.272088;
   double M = 16.014701 * kU - 6 * kMe;
   double Et = Eb + M, pb = std::sqrt(Et * Et - M * M), Ep = ke + kMp, pp = std::sqrt(Ep * Ep - kMp * kMp),
          E4 = Et + kMp - Ep, px = -pp * std::sin(th), pz = pb - pp * std::cos(th), m2 = E4 * E4 - px * px - pz * pz;
   return m2 > 0 ? std::sqrt(m2) - M : -999;
}
} // namespace

void dir3d_score(TString in, int src, double dcaMin = 5.0, double zOff = 1000, int expSense = 0)
{
   TFile f(in);
   auto *t = (TTree *)f.Get("d3");
   int n, ok, spyV, spyD, sen[2];
   double rhoMax, geoTh, tx, ty, tz, tth, ex[2], ey[2], ezd[2], dca[2], R[2], th[2], br[2], qe[2];
   t->SetBranchAddress("n", &n); t->SetBranchAddress("ok", &ok); t->SetBranchAddress("rhoMax", &rhoMax);
   t->SetBranchAddress("geoTh", &geoTh); t->SetBranchAddress("spyV", &spyV); t->SetBranchAddress("spyD", &spyD);
   t->SetBranchAddress("tx", &tx); t->SetBranchAddress("ty", &ty); t->SetBranchAddress("tz", &tz);
   t->SetBranchAddress("tth", &tth); t->SetBranchAddress("ex", ex); t->SetBranchAddress("ey", ey);
   t->SetBranchAddress("ezd", ezd); t->SetBranchAddress("dca", dca); t->SetBranchAddress("R", R);
   t->SetBranchAddress("th", th); t->SetBranchAddress("br", br); t->SetBranchAddress("qe", qe);
   t->SetBranchAddress("sen", sen);

   // ---- pass 1 (sim): which z convention puts an end on the truth vertex ----
   int conv = -1; // z_lab = conv*z_digi + (conv<0 ? zOff : 0)
   if (src == 1) {
      long hit[2] = {0, 0}, tot = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!ok || tz < -900 || rhoMax < 30) continue;
         ++tot;
         for (int c = 0; c < 2; ++c) {
            double s = c ? 1 : -1, o = c ? 0 : zOff, dm = 1e9;
            for (int k = 0; k < 2; ++k)
               dm = std::min(dm, std::sqrt(std::pow(ex[k] - tx, 2) + std::pow(ey[k] - ty, 2) +
                                           std::pow(s * ezd[k] + o - tz, 2)));
            if (dm < 40) ++hit[c];
         }
      }
      conv = hit[1] > hit[0] ? 1 : -1;
      printf("z convention vs truth (end within 40 mm of the generated vertex, %ld tracks):\n"
             "  z_lab = %g - z_digi : %ld   z_lab = z_digi : %ld   -> using %s\n",
             tot, zOff, hit[0], hit[1], conv < 0 ? "z_lab = zOff - z_digi" : "z_lab = z_digi");
   }
   auto zlab = [&](double zd) { return conv < 0 ? zOff - zd : zd; };

   // ---- referee ----
   auto referee = [&](int &vEnd, double &thLab) -> bool {
      if (!ok) return false;
      if (src == 1) {
         if (tz < -900 || rhoMax < 30) return false;
         double d[2];
         for (int k = 0; k < 2; ++k)
            d[k] = std::sqrt(std::pow(ex[k] - tx, 2) + std::pow(ey[k] - ty, 2) + std::pow(zlab(ezd[k]) - tz, 2));
         if (std::min(d[0], d[1]) > 40 || std::max(d[0], d[1]) < 60) return false; // undefined question
         vEnd = d[0] < d[1] ? 0 : 1;
         thLab = tth;
         return true;
      }
      bool on[2];
      double tl[2];
      for (int k = 0; k < 2; ++k) {
         tl[k] = conv < 0 ? 180 - th[k] : th[k];
         double p = 299.792458 * br[k], ke = std::sqrt(p * p + 938.272 * 938.272) - 938.272;
         on[k] = br[k] > 0 && std::abs(exPP(185, tl[k] * TMath::DegToRad(), ke)) < 2 && tl[k] < 90;
      }
      if (on[0] == on[1]) return false;
      vEnd = on[0] ? 0 : 1;
      thLab = tl[vEnd];
      return true;
   };

   // ---- sense calibration: sign seen travelling inward from the TRUE vertex end ----
   if (expSense == 0) {
      long pos = 0, neg = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         int v; double tl;
         if (!referee(v, tl)) continue;
         (sen[v] > 0 ? pos : neg)++;
      }
      expSense = pos >= neg ? 1 : -1;
      printf("rotation sense from the true vertex end: +1 %ld  -1 %ld  -> expected %+d (purity %.1f %%)\n", pos, neg,
             expSense, 100.0 * std::max(pos, neg) / std::max(1L, pos + neg));
   }

   // ---- votes ----
   const int NV = 8;
   const char *vn[NV] = {"rho", "dca", "sense", "radius", "bragg", "spyral", "atpra", "combo"};
   auto votes = [&](int *v) {
      double r0 = std::hypot(ex[0], ey[0]), r1 = std::hypot(ex[1], ey[1]);
      v[0] = r0 == r1 ? -1 : (r0 < r1 ? 0 : 1);
      v[1] = std::abs(dca[0] - dca[1]) < dcaMin ? -1 : (dca[0] < dca[1] ? 0 : 1);
      v[2] = sen[0] == sen[1] ? -1 : (sen[0] == expSense ? 0 : 1);
      v[3] = R[0] == R[1] ? -1 : (R[0] > R[1] ? 0 : 1);
      v[4] = qe[0] == qe[1] ? -1 : (qe[0] < qe[1] ? 0 : 1);
      int lo = ezd[0] < ezd[1] ? 0 : 1;
      v[5] = !spyV ? -1 : (spyD == 0 ? lo : 1 - lo);
      v[6] = std::abs(ezd[0] - ezd[1]) < 1e-9 ? -1 : (geoTh > 90 ? lo : 1 - lo);
      v[7] = v[1] >= 0 ? v[1] : v[2];
   };

   // bins: theta_lab of the referee, and hit count
   const int NB = 6;
   const char *bn[NB] = {"all", "th<80", "80-100", "th>100", "n<150", "n>=500"};
   long good[NB][NV] = {{0}}, dec[NB][NV] = {{0}}, tot[NB] = {0};
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      int v0; double tl;
      if (!referee(v0, tl)) continue;
      int v[NV];
      votes(v);
      bool inb[NB] = {true, tl < 80, tl >= 80 && tl <= 100, tl > 100, n < 150, n >= 500};
      for (int b = 0; b < NB; ++b) {
         if (!inb[b]) continue;
         ++tot[b];
         for (int k = 0; k < NV; ++k)
            if (v[k] >= 0) { ++dec[b][k]; if (v[k] == v0) ++good[b][k]; }
      }
   }
   printf("\naccuracy %% (decided %%) per vote; referee = %s\n", src == 1 ? "MC truth" : "elastic locus");
   printf("%-8s %7s", "bin", "N");
   for (int k = 0; k < NV; ++k) printf(" %14s", vn[k]);
   printf("\n");
   for (int b = 0; b < NB; ++b) {
      printf("%-8s %7ld", bn[b], tot[b]);
      for (int k = 0; k < NV; ++k)
         printf("    %5.1f (%3.0f)", dec[b][k] ? 100.0 * good[b][k] / dec[b][k] : 0, tot[b] ? 100.0 * dec[b][k] / tot[b] : 0);
      printf("\n");
   }
}
