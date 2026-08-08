/// @file apply_theta_corr_C14.C
/// @brief Apply the explorer's empirical theta correction to a cache and recompute Ex / theta_cm.
///
/// The browser explorer's "kinematics correction" block draws theta - slope*(KE - pivot) on the
/// KE vs theta_lab map, but it does NOT propagate: the page's Ex, theta_cm and g.s. fit are still
/// built from the raw theta. This macro makes the correction real, so the corrected angles reach
/// the excitation energy and the angular distribution.
///
///     theta' = theta - (360/N) * (KE - pivot)          [degrees, N is the explorer's kcDenom]
///
/// Both Ex and theta_cm are recomputed from theta' with the same two-body expressions ex_C14.C
/// uses, so the output cache is drop-in: identical column names, and every downstream macro
/// (elastic_sideband_C14.C, the acceptance, the FRESCO comparison) works on it unchanged.
///
/// ACCEPTANCE NOTE. The simulated acceptance is a function of TRUE theta_cm and was built from a
/// sim that does NOT carry the data's +8 % KE bias. The correction moves the DATA's theta_cm
/// toward truth, so pairing corrected data with the uncorrected acceptance is more consistent
/// than the uncorrected pairing was -- not less. The correction must NOT be applied to the sim.
///
///   root -b -q 'apply_theta_corr_C14.C("plots/proton_kin_300gfx.root","_300gfx_tc",1125,3.5)'

#include <tuple>

static double tc_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// verbatim from ex_C14.C: returns {Ex, theta_cm [deg]}
static std::tuple<double, double> tc_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * tc_om2(s, m1 * m1, m2 * m2) * tc_om2(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg <= 0)
      return {NAN, NAN};
   double m4_ex = std::sqrt(arg);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double ct = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
               (tc_om2(s, m1 * m1, m2 * m2) * tc_om2(s, m3 * m3, m4_ex * m4_ex));
   ct = std::max(-1.0, std::min(1.0, ct));
   return {Ex, (TMath::Pi() - std::acos(ct)) * TMath::RadToDeg()};
}

void apply_theta_corr_C14(TString inCache = "plots/proton_kin_300gfx.root", TString outTag = "_300gfx_tc",
                          Double_t kcDenom = 1125.0, Double_t kcPivot = 3.5, Double_t Ebeam = 161.0,
                          Double_t mEjectAmu = 1.007825, Double_t mResidAmu = 14.003242)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401;
   const double m1 = 14.003242 * u, m2 = 1.007825 * u, m3 = mEjectAmu * u, m4 = mResidAmu * u;
   const double slopeDeg = 360.0 / kcDenom; // the explorer stores N in 2*pi/N rad/MeV

   TFile *fi = TFile::Open(here + "/" + inCache);
   if (!fi || fi->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", inCache.Data());
      return;
   }
   TTree *t = (TTree *)fi->Get("pk");
   if (!t) {
      printf("\033[1;31mno tree `pk`\033[0m\n");
      return;
   }
   float ke, th, vz = 0, thcm, ex, c2n;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("thcm", &thcm);
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("chi2ndf", &c2n);
   const bool hasVz = t->GetBranch("vertexz") != nullptr;
   if (hasVz)
      t->SetBranchAddress("vertexz", &vz);

   TString outPath = here + "/plots/proton_kin" + outTag + ".root";
   TFile fo(outPath, "RECREATE");
   auto *nt = new TNtuple("pk", "proton kinematics (theta-corrected)", "ke:theta:vertexz:thcm:ex:chi2ndf");

   printf("\n===== theta correction: theta' = theta - %.4f deg/MeV * (KE - %.2f MeV) =====\n", slopeDeg, kcPivot);
   printf("      (explorer N = %.1f in 2pi/N rad/MeV, pivot = %.2f MeV, Ebeam = %.1f)\n", kcDenom, kcPivot, Ebeam);

   Long64_t N = t->GetEntries();
   long kept = 0, dropped = 0;
   double sumShift = 0, maxShift = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      double dth = slopeDeg * (ke - kcPivot);
      double thC = th - dth;
      if (!(thC > 0.05 && thC < 179.95)) { // a correction that walks the angle out of range
         ++dropped;
         continue;
      }
      auto [exC, cmC] = tc_kine(m1, m2, m3, m4, Ebeam, thC * TMath::DegToRad(), ke);
      if (!std::isfinite(exC) || !std::isfinite(cmC)) {
         ++dropped;
         continue;
      }
      sumShift += std::fabs(dth);
      maxShift = std::max(maxShift, std::fabs(dth));
      nt->Fill(ke, thC, vz, cmC, exC, c2n);
      ++kept;
   }
   printf("  %lld tracks in -> %ld written, %ld dropped (corrected angle unphysical)\n", N, kept, dropped);
   printf("  mean |dtheta| %.2f deg, max %.2f deg\n", kept ? sumShift / kept : 0, maxShift);
   nt->Write();
   fo.Close();
   fi->Close();
   printf("  wrote %s\n\n", outPath.Data());
}
