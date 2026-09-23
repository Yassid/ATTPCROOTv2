// exsel_dp.C -- reproduce the BROWSER EXPLORER's (d,p) selection and Ex offline.
//
// WHY THIS EXISTS: the explorer page carries the adopted selection, but it is a static HTML
// page with no persistence -- the only durable record is the params JSON it exports. This
// macro is the headless twin of that page, so a spectrum can be fitted instead of eyeballed.
//
// EVERY formula below is COPIED from the generated page, not re-derived:
//   kine2b/omega2 ............ the page's two-body invariant-mass kinematics
//   thCorr ................... th - (360/kcDenom)*(ke + keOff - kcPivot)
//   pass() ................... c2 <= chi2 && thCorr in [thLo,thHi] && ke+keOff in [keLo,keHi]
//                              && vz in [vzLo,vzHi]
//   keOff .................... added to ke BEFORE the cut, the theta correction and the
//                              kinematics alike, exactly as add_keoff.py patches it
// Re-deriving the Ex kinematics by hand is a 0.511 MeV/electron mass-convention trap
// (nuclear light / atomic heavy); do not "tidy" the mass constants.
//
// Defaults are the GOOD_explorer_params_16C_d_p_17C_..._2026-09-07-23-36-03.json values.
//
//   root -l 'exsel_dp.C("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root")'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <cstdio>
#include <vector>

namespace exsel {

// masses: the page's CFG, amu -> MeV with the page's own U
constexpr double U = 931.49401;
constexpr double M1 = 16.0147013 * U;      // 16C beam      (atomic)
constexpr double M2 = 2.0135532 * U;       // d target      (nuclear)
constexpr double M3 = 1.00727646688 * U;   // p ejectile    (nuclear)
constexpr double M4 = 17.0225864 * U;      // 17C residual  (atomic)

inline double omega2(double x, double y, double z)
{
   return TMath::Sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics -> (Ex, theta_cm in deg); Ex = NaN below threshold
inline bool kine2b(double Kp, double thl_rad, double Ke, double &ex, double &tcm_deg)
{
   const double Et1 = Kp + M1, Et3 = Ke + M3, Et4 = Et1 + M2 - Et3;
   const double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   const double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   const double a = (TMath::Cos(thl_rad) * omega2(s, M1 * M1, M2 * M2) * omega2(uu, M2 * M2, M3 * M3) -
                     (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                       (2 * M2 * M2) +
                    s + uu - M2 * M2;
   if (a < 0) return false;
   const double m4x = TMath::Sqrt(a);
   ex = m4x - M4;
   const double t = M2 * M2 + m4x * m4x - 2 * M2 * Et4;
   double arg = (s * s + s * (2 * t - M1 * M1 - M2 * M2 - M3 * M3 - m4x * m4x) +
                 (M1 * M1 - M2 * M2) * (M3 * M3 - m4x * m4x)) /
                (omega2(s, M1 * M1, M2 * M2) * omega2(s, M3 * M3, m4x * m4x));
   arg = TMath::Max(-1.0, TMath::Min(1.0, arg));
   tcm_deg = (TMath::Pi() - TMath::ACos(arg)) * 180.0 / TMath::Pi();
   return true;
}

struct Cuts {
   double ebeam = 184.25;
   double chi2 = 13.5;
   double thLo = 89, thHi = 161;
   double keLo = 1.2, keHi = 1000;
   double vzLo = -100, vzHi = 940;
   double keOff = 0.25;
   bool kcOn = true;
   double kcDenom = 2250.0, kcPivot = 3.0;
   // Ex histogram, as the page has it
   int exBins = 90;
   double exLo = -5, exHi = 9;
};

inline double kcSlopeDeg(const Cuts &s)
{
   return (TMath::Finite(s.kcDenom) && TMath::Abs(s.kcDenom) > 1e-9) ? 360.0 / s.kcDenom : 0.0;
}
inline double thCorr(const Cuts &s, double th, double ke)
{
   return s.kcOn ? th - kcSlopeDeg(s) * (ke - s.kcPivot) : th;
}

} // namespace exsel

void exsel_dp(const char *fname = "/mnt/f/a1975/caches/.explorer/exp_dp_find805.root",
              const char *outname = "", double chi2 = 13.5, double thLo = 89, double thHi = 161,
              double keLo = 1.2, double keOff = 0.25, double ebeam = 184.25)
{
   using namespace exsel;
   Cuts s;
   s.chi2 = chi2; s.thLo = thLo; s.thHi = thHi; s.keLo = keLo; s.keOff = keOff; s.ebeam = ebeam;

   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", fname); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("ERROR: no tree 'pk' in %s\n", fname); return; }

   float ke = 0, th = 0, vz = 0, c2 = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);

   TH1D *hEx = new TH1D("hEx", "16C(d,p)17C  E_{x};E_{x} [MeV];counts / bin", s.exBins, s.exLo, s.exHi);
   TH2D *hExTh = new TH2D("hExTh", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]",
                          60, 0, 180, s.exBins, s.exLo, s.exHi);
   hEx->Sumw2();

   const Long64_t N = t->GetEntries();
   Long64_t nPass = 0, nKine = 0, nInHist = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      const double k = ke + s.keOff;
      const double thc = thCorr(s, th, k);
      if (!(c2 <= s.chi2)) continue;
      if (!(thc >= s.thLo && thc <= s.thHi)) continue;
      if (!(k >= s.keLo && k <= s.keHi)) continue;
      if (!(vz >= s.vzLo && vz <= s.vzHi)) continue;
      ++nPass;
      double ex = 0, tcm = 0;
      if (!kine2b(s.ebeam, thc * TMath::Pi() / 180.0, k, ex, tcm)) continue;
      ++nKine;
      hEx->Fill(ex);
      hExTh->Fill(tcm, ex);
      if (ex >= s.exLo && ex <= s.exHi) ++nInHist;
   }

   printf("\n=== exsel_dp: %s ===\n", fname);
   printf("  Ebeam %.2f   chi2ndf <= %.2f   thCorr [%.0f,%.0f]   ke+%.2f in [%.2f,%.0f]   vz [%.0f,%.0f]\n",
          s.ebeam, s.chi2, s.thLo, s.thHi, s.keOff, s.keLo, s.keHi, s.vzLo, s.vzHi);
   printf("  theta correction %s: slope %.4f deg/MeV about pivot %.2f MeV\n",
          s.kcOn ? "ON" : "OFF", kcSlopeDeg(s), s.kcPivot);
   printf("  tracks in ntuple      %lld\n", N);
   printf("  pass the cuts         %lld  (%.1f%%)\n", nPass, 100.0 * nPass / N);
   printf("  kinematics solvable   %lld\n", nKine);
   printf("  inside Ex [%.0f,%.0f]   %lld\n", s.exLo, s.exHi, nInHist);
   printf("  bin width             %.4f MeV\n", (s.exHi - s.exLo) / s.exBins);

   // where the counts sit, in the regions that matter for 17C
   const double Sn = 0.729;
   auto range = [&](double a, double b) {
      return hEx->Integral(hEx->FindBin(a + 1e-6), hEx->FindBin(b - 1e-6));
   };
   printf("\n  BOUND      Ex < %.3f (Sn)        %8.0f\n", Sn, range(s.exLo, Sn));
   printf("    of which Ex < -5 (unphysical) %8.0f  [underflow]\n", hEx->GetBinContent(0));
   printf("  UNBOUND    %.3f < Ex < 9        %8.0f\n", Sn, range(Sn, 9));
   printf("    2.0 - 3.5                     %8.0f\n", range(2.0, 3.5));
   printf("    3.5 - 5.0                     %8.0f\n", range(3.5, 5.0));
   printf("    5.0 - 7.0                     %8.0f\n", range(5.0, 7.0));
   printf("    7.0 - 9.0                     %8.0f\n", range(7.0, 9.0));

   if (outname && outname[0]) {
      TFile *fo = TFile::Open(outname, "RECREATE");
      hEx->Write();
      hExTh->Write();
      fo->Close();
      printf("\n  wrote %s\n", outname);
   }
   printf("\n");
}
