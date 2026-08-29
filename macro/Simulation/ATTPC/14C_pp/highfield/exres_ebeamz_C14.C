/// @file exres_ebeamz_C14.C
/// @brief What the Ex resolution would be if the beam energy at the vertex were used instead of
/// one constant for the whole chamber.
///
/// WHY. The campaign reproduces the adopted a1954 analysis, which reconstructs every event with
/// the same beam energy (159.75 MeV) whatever the vertex. The 14C beam loses about 12 MeV
/// crossing 1 m of H2 at 300 torr and the vertex is uniform in z, so that single number is wrong
/// by up to several MeV in a way that is perfectly correlated with a measured quantity. Rebuild
/// Ex from PERFECT truth kinematics under that assumption and the spread is already 0.25 MeV --
/// the same as the fully reconstructed spread. In other words the measurement as analysed is not
/// detector-limited at all, and no pad pitch or field can show up until this is dealt with.
///
/// WHAT THIS MACRO DOES, in three passes over the residual tree written by ex_res_C14_hf.C:
///   1. CALIBRATE. For each event the true (KE, theta) of the recoil proton are known and the
///      true excitation is known, so there is exactly one beam energy that reproduces it. Solve
///      for it by bisection, event by event, and fit the result against the true vertex z. That
///      curve IS the beam energy profile of this simulation, extracted from it rather than
///      assumed, and its scatter about the fit is the straggling -- the part that stays
///      irreducible however well the vertex is measured.
///   2. IDEAL. Rebuild Ex with the fitted profile evaluated at the TRUE vertex z. This is the
///      floor: perfect vertex, deterministic energy loss, only straggling and tracking left.
///   3. REALISTIC. Rebuild Ex with the profile evaluated at the RECONSTRUCTED vertex z, which is
///      what an experiment can actually do.
///
/// The three numbers together say how much of the gap is the method and how much is the detector.
///
///   root -b -q 'exres_ebeamz_C14.C("/mnt/f/a1954_C14_hf/b400_2mm/exres_gs_s7101_b400_2mm.root",0.0)'

#include <algorithm>
#include <vector>

static double ez_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// Ex from the ejectile (KE, theta_lab) for a given beam energy -- the same expression the
/// analysis uses, so nothing here is a different kinematics convention
static double ez_ex(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * ez_omega2(s, m1 * m1, m2 * m2) * ez_omega2(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg <= 0)
      return NAN;
   return std::sqrt(arg) - m4;
}

static double ez_q(std::vector<double> v, double p)
{
   if (v.empty())
      return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static void ez_report(const char *label, std::vector<double> v)
{
   if (v.size() < 20) {
      printf("  %-34s %6zu       --        --        --\n", label, v.size());
      return;
   }
   double a = ez_q(v, .25), b = ez_q(v, .50), c = ez_q(v, .75);
   printf("  %-34s %6zu  %+8.3f  %8.3f  %8.3f\n", label, v.size(), b, c - a, (c - a) / 1.349);
}

void exres_ebeamz_C14(TString resFile, Double_t resEx = 0.0, Double_t Ebeam0 = 159.75, Double_t thMin = 20.,
                      Double_t thMax = 90.)
{
   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;
   const double m_resid = m_C14 + resEx;

   TFile *fr = TFile::Open(resFile);
   if (!fr || fr->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", resFile.Data());
      return;
   }
   TTree *t = (TTree *)fr->Get("res");
   if (!t) {
      printf("\033[1;31mno res tree\033[0m\n");
      return;
   }
   double exReco, exTrue, thTrue, thReco, keTrue, keReco, zTrue, zReco;
   t->SetBranchAddress("exReco", &exReco);
   t->SetBranchAddress("exTrue", &exTrue);
   t->SetBranchAddress("thTrue", &thTrue);
   t->SetBranchAddress("thReco", &thReco);
   t->SetBranchAddress("keTrue", &keTrue);
   t->SetBranchAddress("keReco", &keReco);
   t->SetBranchAddress("zTrue", &zTrue);
   t->SetBranchAddress("zReco", &zReco);

   // ---- pass 1 : the beam energy each truth event requires ---------------------------------
   auto *pEb = new TProfile("pEb", ";z_{true} [mm];E_{beam} required [MeV]", 25, 0, 1000, "s");
   std::vector<double> ebv, zv;
   Long64_t N = t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (thTrue < thMin || thTrue >= thMax)
         continue;
      // bisection on the beam energy: Ex(E) is monotonic over any sane bracket
      double lo = 100., hi = 200.;
      double flo = ez_ex(m_C14, m_p, m_p, m_resid, lo, thTrue * TMath::DegToRad(), keTrue);
      double fhi = ez_ex(m_C14, m_p, m_p, m_resid, hi, thTrue * TMath::DegToRad(), keTrue);
      if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0)
         continue;
      for (int it = 0; it < 60; ++it) {
         double mid = 0.5 * (lo + hi);
         double f = ez_ex(m_C14, m_p, m_p, m_resid, mid, thTrue * TMath::DegToRad(), keTrue);
         if (std::isnan(f))
            break;
         if (f * flo <= 0) {
            hi = mid;
            fhi = f;
         } else {
            lo = mid;
            flo = f;
         }
      }
      double eb = 0.5 * (lo + hi);
      if (eb < 105 || eb > 195)
         continue;
      ebv.push_back(eb);
      zv.push_back(zTrue);
      pEb->Fill(zTrue, eb);
   }
   if (ebv.size() < 100) {
      printf("\033[1;31monly %zu events solved -- nothing to fit\033[0m\n", ebv.size());
      return;
   }
   // quadratic in z: the loss per unit length rises as the beam slows, so a straight line leaves
   // a visible curvature behind
   TGraph g((int)ebv.size(), zv.data(), ebv.data());
   TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
   fEb.SetParameters(Ebeam0, -0.012, 0.);
   g.Fit(&fEb, "QN");
   printf("\n=== %s ===\n", gSystem->BaseName(resFile.Data()));
   printf("beam-energy profile extracted from truth:  E(z) = %.3f %+.5f z %+.3e z^2  MeV  (z in mm)\n",
          fEb.GetParameter(0), fEb.GetParameter(1), fEb.GetParameter(2));
   printf("   E(0) = %.2f   E(500) = %.2f   E(1000) = %.2f MeV   -> %.2f MeV lost across the drift\n",
          fEb.Eval(0), fEb.Eval(500), fEb.Eval(1000), fEb.Eval(0) - fEb.Eval(1000));
   std::vector<double> strag;
   for (size_t i = 0; i < ebv.size(); ++i)
      strag.push_back(ebv[i] - fEb.Eval(zv[i]));
   printf("   straggling about the profile: IQR %.3f MeV (sigma %.3f)\n", ez_q(strag, .75) - ez_q(strag, .25),
          (ez_q(strag, .75) - ez_q(strag, .25)) / 1.349);

   // ---- passes 2 and 3 : rebuild Ex --------------------------------------------------------
   std::vector<double> dConst, dTruthOnly, dIdeal, dReal, dRealTruthKin;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (thTrue < thMin || thTrue >= thMax)
         continue;
      dConst.push_back(exReco - resEx);
      dTruthOnly.push_back(exTrue);
      double ebT = fEb.Eval(zTrue), ebR = fEb.Eval(zReco);
      double a = ez_ex(m_C14, m_p, m_p, m_C14, ebT, thReco * TMath::DegToRad(), keReco);
      double b = ez_ex(m_C14, m_p, m_p, m_C14, ebR, thReco * TMath::DegToRad(), keReco);
      double c = ez_ex(m_C14, m_p, m_p, m_C14, ebR, thTrue * TMath::DegToRad(), keTrue);
      if (!std::isnan(a))
         dIdeal.push_back(a - resEx);
      if (!std::isnan(b))
         dReal.push_back(b - resEx);
      if (!std::isnan(c))
         dRealTruthKin.push_back(c - resEx);
   }

   printf("\n  %-34s %6s  %8s  %8s  %8s\n", "Ex reconstruction", "n", "median", "IQR", "sigma");
   ez_report("constant E_beam (as analysed)", dConst);
   ez_report("  same, but perfect tracking", dTruthOnly);
   ez_report("E_beam(z_true), reco kinematics", dIdeal);
   ez_report("E_beam(z_reco), reco kinematics", dReal);
   ez_report("  same, but perfect tracking", dRealTruthKin);
   printf("\nThe last row is the METHOD floor at a reconstructed vertex: whatever it says, no pad\n"
          "pitch or field can go under it. The gap between it and the row above is the detector.\n");
   printf("\nebeamz done\n\n");
   fr->Close();
}
