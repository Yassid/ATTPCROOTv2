/// @file dp_summary_C17.C
/// @brief The 17C(d,p)18C result: excitation-energy resolution over the FULL lab range, the
///        resolution floor beside it, the vertex beam-energy correction, and whether the 18C
///        levels separate.
///
///   root -b -q 'dp_summary_C17.C'
///   root -b -q 'dp_summary_C17.C("/mnt/f/C17dp/b285_attpc/")'
///
/// WHY THIS EXISTS RATHER THAN THE PER-SAMPLE TABLE. ex_res_C14_hf.C hard-codes theta_lab slices
/// 20-90 deg: it was written for 14C(p,p'), where the recoil proton obeys
/// theta_lab = (180-theta_cm)/2 and 90 deg IS the physical limit. A (d,p) proton has no such
/// limit, and in THIS channel the transfer peak (theta_cm 2-40) sits at theta_lab 100-175 deg --
/// entirely outside those slices. Its per-slice table is therefore blind exactly where the physics
/// is. (Its "ALL" row is complete, and its flat "res" TTree holds every accepted event, which is
/// what this macro re-bins. Nothing has to be re-run.)
///
/// WHAT IT REPORTS, and why each one is here rather than a single number:
///
///  1. sigma(Ex) per theta_lab slice across the WHOLE range, as IQR/1.349. A robust width, not an
///     RMS and not a fit: fit-failure tails are heavy enough in this pipeline to fake large
///     effects, and a median-based width does not care.
///  2. THE RESOLUTION FLOOR beside every measured width -- the same inversion applied to the TRUTH
///     (keTrue, thTrue) at the same constant beam energy. The floor is what the METHOD costs when
///     the detector is perfect, so floor ~ measured means the detector is invisible and no
///     detector change can help, while floor << measured means there is something to chase.
///     Quoting a resolution without it cannot distinguish those two cases.
///  3. THE TAIL FRACTION (|Ex - Ex_true| > 1 MeV). A width alone cannot tell a clean measurement
///     from one with a bad systematic sitting under a narrow core.
///  4. THE VERTEX BEAM-ENERGY CORRECTION. The beam loses 14.8 MeV crossing the chamber, so a
///     constant Ebeam is wrong by up to +-7 MeV depending on where the reaction happened. In
///     14C(d,p) taking Ebeam at the reconstructed vertex instead took backward sigma(Ex) from
///     0.178 to 0.064 MeV -- it is the single largest term and it costs nothing but software.
///     E_beam(z) is MEASURED here from the truth, not assumed: for each event the beam energy that
///     makes the known Ex close exactly is solved for, then fitted against z.
///  5. LEVEL SEPARATION, as Delta(Ex) / (sigma_1 + sigma_2) for each adjacent pair of 18C bound
///     states. The tightest is 1.588 -> 2.515, a 927 keV gap, against the ~300 keV the proposal
///     quotes for the AT-TPC.
#include <algorithm>
#include <map>
#include <vector>

static double s17_omega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

/// Two-body inversion, identical in form to acceptance_C14.C:acc_kine -- given the measured
/// ejectile (theta_lab, KE) and a beam energy, return the residual excitation energy.
static double s17_ex(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * s17_omega2(s, m1 * m1, m2 * m2) * s17_omega2(u, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                         (2 * m2 * m2) +
                      s + u - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

/// robust width: the interquartile range over 1.349, which equals sigma for a gaussian but is
/// immune to the fit-failure tail. Returns -1 on too little data.
static double s17_iqrSigma(std::vector<double> v)
{
   if (v.size() < 20)
      return -1;
   std::sort(v.begin(), v.end());
   const double q1 = v[(size_t)(0.25 * v.size())], q3 = v[(size_t)(0.75 * v.size())];
   return (q3 - q1) / 1.349;
}
static double s17_median(std::vector<double> v)
{
   if (v.empty())
      return 0;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}

void dp_summary_C17(TString dir = "/mnt/f/C17dp/b285_attpc/", Double_t EbeamConst = 135.0, Double_t chi2Cut = 5.0,
                    Double_t mBeamAmu = 17.0225787, Double_t mTgtAmu = 2.0141018, Double_t mEjAmu = 1.007825,
                    Double_t mResAmu = 18.0267519)
{
   const double u = 931.49401;
   const double m1 = mBeamAmu * u, m2 = mTgtAmu * u, m3 = mEjAmu * u, m4 = mResAmu * u;

   const int nL = 4;
   const double ExGen[nL] = {0.0, 1.588, 2.515, 3.972};
   const char *stTag[nL] = {"gs", "ex1588", "ex2515", "ex3972"};
   const char *stName[nL] = {"0+ g.s.", "2+ 1.588", "(2+) 2.515", "(2,3)+ 3.972"};

   // full lab range, in 15 deg slices -- the point of this macro
   const int nSl = 12;
   const double slLo[nSl] = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165};
   const double slHi[nSl] = {15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180};

   printf("\n\033[1;33m########## 17C(d,p)18C -- D2 300 torr, B = 2.85 T, AT-TPC pad plane ##########\033[0m\n");
   printf("  constant-Ebeam analysis at %.2f MeV, chi2/ndf < %.1f\n", EbeamConst, chi2Cut);

   // per level: the residual samples, constant-Ebeam and vertex-corrected
   std::vector<double> resC[nL], resV[nL], floorC[nL];
   std::vector<double> resCsl[nL][nSl], floorCsl[nL][nSl], resVsl[nL][nSl];
   std::vector<double> resPeakC[nL], resPeakV[nL], floorPeakC[nL];
   long nTot[nL] = {0}, nTailC[nL] = {0}, nTailV[nL] = {0};
   double ebz_a = 0, ebz_b = 0; // E_beam(z) = a + b*z, measured from truth
   bool haveEbz = false;

   for (int l = 0; l < nL; ++l) {
      // find the sample file for this level, whatever seed it carries
      TString pat = TString(dir) + "exres_" + stTag[l] + "_s*_b285_attpc.root";
      TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-12s MISSING (%s)\033[0m\n", stTag[l], pat.Data());
         continue;
      }
      TFile *f = TFile::Open(found);
      if (!f || f->IsZombie()) {
         printf("\033[1;31m  cannot open %s\033[0m\n", found.Data());
         continue;
      }
      TTree *t = (TTree *)f->Get("res");
      if (!t) {
         printf("\033[1;31m  no res tree in %s\033[0m\n", found.Data());
         continue;
      }
      double exR, exT, thT, thR, keT, keR, cmT, c2n, zT, zR;
      t->SetBranchAddress("exReco", &exR);
      t->SetBranchAddress("exTrue", &exT);
      t->SetBranchAddress("thTrue", &thT);
      t->SetBranchAddress("thReco", &thR);
      t->SetBranchAddress("keTrue", &keT);
      t->SetBranchAddress("keReco", &keR);
      t->SetBranchAddress("cmTrue", &cmT);
      t->SetBranchAddress("chi2ndf", &c2n);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);

      // ---- FIRST PASS on the ground state: MEASURE E_beam(z) from the truth. For each event the
      // ---- beam energy that makes the GENERATED Ex close exactly is solved for by bisection,
      // ---- then fitted against the true vertex z. Nothing is assumed from a stopping-power table.
      //
      // NOTE ON THE REFERENCE, which is easy to get wrong and silently self-fulfilling: the
      // reference here is ExGen[l], the excitation the generator was ASKED for -- NOT the tree's
      // exTrue. ex_res_C14_hf.C:162 builds exTrue by putting the TRUTH kinematics through this
      // same inversion at this same constant beam energy (with the excited residual mass, so it
      // sits near 0). Solving for the beam energy that reproduces exTrue would therefore return
      // EbeamConst exactly, by construction, and a "floor" measured against exTrue would be
      // identically zero. Both happened on the first version of this macro.
      if (!haveEbz && l == 0) {
         TGraph g;
         const double tgt = ExGen[l];
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = EbeamConst - 40, hi = EbeamConst + 40;
            if ((s17_ex(m1, m2, m3, m4, lo, thT, keT) - tgt) * (s17_ex(m1, m2, m3, m4, hi, thT, keT) - tgt) > 0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((s17_ex(m1, m2, m3, m4, lo, thT, keT) - tgt) * (s17_ex(m1, m2, m3, m4, mid, thT, keT) - tgt) <= 0)
                  hi = mid;
               else
                  lo = mid;
            }
            g.SetPoint(g.GetN(), zT, 0.5 * (lo + hi));
         }
         if (g.GetN() > 100) {
            TF1 lin("lin", "pol1");
            g.Fit(&lin, "QN");
            ebz_a = lin.GetParameter(0);
            ebz_b = lin.GetParameter(1);
            haveEbz = true;
            printf("\n\033[1;33m===== beam energy at the vertex, MEASURED from truth =====\033[0m\n");
            printf("  E_beam(z) = %.4f %+.6f * z[mm]   (n = %d)\n", ebz_a, ebz_b, g.GetN());
            printf("  -> %.2f MeV at z = 0, %.2f MeV at z = 1000 mm; loss over the metre %.2f MeV\n", ebz_a,
                   ebz_a + 1000 * ebz_b, -1000 * ebz_b);
            printf("  (CATIMA predicts 142.29 -> 127.49, i.e. 14.80 MeV; the constant used is %.2f)\n", EbeamConst);
         }
      }

      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (c2n >= chi2Cut)
            continue;
         ++nTot[l];
         // THE RESIDUAL, against the GENERATED excitation -- the same quantity
         // ex_res_C14_hf.C:208 computes as `exR - resEx`. exR is built with the GROUND-STATE
         // residual mass, because an experiment does not know the level in advance, so it should
         // come out at ExGen and this difference is the error.
         const double rC = exR - ExGen[l];
         // THE FLOOR, free: exTrue is already the truth kinematics through this same inversion at
         // this same constant beam energy, using the EXCITED residual mass, so it sits at 0 for a
         // perfect method and its offset from 0 IS the method term. See the note above.
         const double fC = exT;
         // vertex-corrected: the measured kinematics, but the beam energy at the RECONSTRUCTED z
         double rV = rC;
         if (haveEbz)
            rV = s17_ex(m1, m2, m3, m4, ebz_a + ebz_b * zR, thR, keR) - ExGen[l];

         resC[l].push_back(rC);
         resV[l].push_back(rV);
         floorC[l].push_back(fC);
         if (std::fabs(rC) > 1.0)
            ++nTailC[l];
         if (std::fabs(rV) > 1.0)
            ++nTailV[l];
         for (int s = 0; s < nSl; ++s)
            if (thT >= slLo[s] && thT < slHi[s]) {
               resCsl[l][s].push_back(rC);
               resVsl[l][s].push_back(rV);
               floorCsl[l][s].push_back(fC);
            }
         // the transfer peak: small theta_cm, where a stripping angular distribution has its yield
         if (cmT >= 2.0 && cmT <= 40.0) {
            resPeakC[l].push_back(rC);
            resPeakV[l].push_back(rV);
            floorPeakC[l].push_back(fC);
         }
      }
   }

   // ---- 1+2+3. resolution vs theta_lab, with the floor beside it -------------------------------
   printf("\n\033[1;33m===== sigma(E_x) vs theta_lab, FULL RANGE (IQR/1.349, MeV) =====\033[0m\n");
   printf("  'floor' = the same inversion on TRUTH kinematics: what the constant-Ebeam METHOD costs\n");
   printf("  with a perfect detector. floor ~ measured means the detector is not the limit.\n\n");
   printf("  %-9s", "theta_lab");
   for (int l = 0; l < nL; ++l)
      printf(" | %18s", stName[l]);
   printf("\n  %-9s", "[deg]");
   for (int l = 0; l < nL; ++l)
      printf(" | %6s %6s %4s", "meas", "floor", "n");
   printf("\n");
   for (int s = 0; s < nSl; ++s) {
      printf("  %3.0f-%-5.0f", slLo[s], slHi[s]);
      for (int l = 0; l < nL; ++l) {
         const double sm = s17_iqrSigma(resCsl[l][s]), sf = s17_iqrSigma(floorCsl[l][s]);
         if (sm < 0)
            printf(" | %6s %6s %4d", "-", "-", (int)resCsl[l][s].size());
         else
            printf(" | %6.3f %6.3f %4d", sm, sf, (int)resCsl[l][s].size());
      }
      printf("\n");
   }

   // ---- 4. the transfer peak, and what the vertex correction buys -------------------------------
   printf("\n\033[1;33m===== the TRANSFER PEAK (theta_cm 2-40 deg, i.e. theta_lab ~100-175) =====\033[0m\n");
   printf("  %-14s %6s %9s %9s %9s %9s %9s\n", "level", "n", "sig_const", "sig_vtx", "floor", "med_const",
          "tail_vtx");
   for (int l = 0; l < nL; ++l) {
      if (resPeakC[l].empty()) {
         printf("  %-14s %6s\n", stName[l], "-");
         continue;
      }
      long tail = 0;
      for (double r : resPeakV[l])
         if (std::fabs(r) > 1.0)
            ++tail;
      printf("  %-14s %6d %9.3f %9.3f %9.3f %9.3f %8.1f %%\n", stName[l], (int)resPeakC[l].size(),
             s17_iqrSigma(resPeakC[l]), s17_iqrSigma(resPeakV[l]), s17_iqrSigma(floorPeakC[l]),
             s17_median(resPeakC[l]), 100.0 * tail / resPeakV[l].size());
   }

   printf("\n\033[1;33m===== ALL accepted tracks =====\033[0m\n");
   printf("  %-14s %7s %9s %9s %9s %9s %9s\n", "level", "n", "sig_const", "sig_vtx", "floor", "tail_const",
          "tail_vtx");
   for (int l = 0; l < nL; ++l) {
      if (!nTot[l]) {
         printf("  %-14s %7s\n", stName[l], "-");
         continue;
      }
      printf("  %-14s %7ld %9.3f %9.3f %9.3f %8.1f %% %8.1f %%\n", stName[l], nTot[l], s17_iqrSigma(resC[l]),
             s17_iqrSigma(resV[l]), s17_iqrSigma(floorC[l]), 100.0 * nTailC[l] / nTot[l],
             100.0 * nTailV[l] / nTot[l]);
   }

   // ---- 5. do the levels separate? ---------------------------------------------------------------
   printf("\n\033[1;33m===== level separation: Delta(Ex) / (sigma_1 + sigma_2) =====\033[0m\n");
   printf("  A pair is comfortably resolved at >~ 2, marginal near 1, blended below.\n");
   printf("  %-28s %9s %11s %11s\n", "pair", "Delta", "peak(vtx)", "all(vtx)");
   for (int l = 0; l + 1 < nL; ++l) {
      const double d = ExGen[l + 1] - ExGen[l];
      const double sp1 = s17_iqrSigma(resPeakV[l]), sp2 = s17_iqrSigma(resPeakV[l + 1]);
      const double sa1 = s17_iqrSigma(resV[l]), sa2 = s17_iqrSigma(resV[l + 1]);
      TString pr = TString::Format("%s -> %s", stName[l], stName[l + 1]);
      printf("  %-28s %9.3f", pr.Data(), d);
      if (sp1 > 0 && sp2 > 0)
         printf(" %11.2f", d / (sp1 + sp2));
      else
         printf(" %11s", "-");
      if (sa1 > 0 && sa2 > 0)
         printf(" %11.2f", d / (sa1 + sa2));
      else
         printf(" %11s", "-");
      printf("\n");
   }
   printf("\n  The tightest pair in 18C is 1.588 -> 2.515, a 927 keV gap; the proposal quotes\n"
          "  ~300 keV as the AT-TPC excitation-energy resolution.\n\n");
}
