/// @file ab_catima_cmp.C
/// @brief Compare two genfit productions that differ only in the CATIMA material backend.
///
/// The observable is KE_fit / KE_Brho, where KE_Brho comes from the Spyral arc fit and is
/// INDEPENDENT of genfit: it is curvature, not a Kalman filter, so it does not move when the
/// material model changes. The ratio therefore measures the fitter against a fixed ruler.
///
/// The comparison is PAIRED. Tracks are matched by (entry, trackID) across the two files and
/// only tracks present in BOTH contribute to the ratio statistics, so a difference cannot be
/// manufactured by the two arms keeping different subsets. Arm-only tracks are counted and
/// reported separately -- with matFallback off, a throw is a real difference between the
/// models and is worth seeing rather than hiding.
///
/// Medians and 16-84 percentiles throughout: the fit-failure tail on this data is heavy enough
/// that a mean and an RMS have faked an entire result before.
///
///   root -b -q 'ab_catima_cmp.C("/mnt/f/a1975/gf_dt_ab_off/run_0031_multifit_genfitter_t.root",
///                               "/mnt/f/a1975/gf_dt_ab_on/run_0031_multifit_genfitter_t.root")'

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace {

constexpr double kMassTriton = 2808.921; // MeV
constexpr double kSplit = 3.5;           // MeV, beta*gamma = 0.05 for a triton

struct Rec {
   double keFit;  // raw fit at the first measurement point
   double keXtr;  // back-extrapolated to the vertex
   double keBrho; // independent: from the Spyral arc fit
   double chi2n;
   double ndf;
   double chi2;
   bool matUsed;
   bool fellBack;
   // ndf < 0 with chi2 == 0 is a COLLAPSED fit: the filter kept essentially no measurement, so
   // ndf = (measurement dims) - 5 went negative. Such a track still carries kinematics and still
   // reaches the ntuple, which is exactly how a fit-failure tail gets mistaken for a resolution.
   bool good() const { return ndf > 0 && chi2 > 0; }
};

double keFromBrho(double brho)
{
   const double p = 299.792458 * brho; // Z = 1, brho in T*m -> MeV/c
   return std::sqrt(p * p + kMassTriton * kMassTriton) - kMassTriton;
}

// key = entry * 10000 + trackID
std::map<long, Rec> harvest(TString file, long &nEvt, long &nFit)
{
   std::map<long, Rec> out;
   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", file.Data());
      return out;
   }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   t->SetBranchAddress("AtPIDEvent", &pe);
   nEvt = t->GetEntries();
   for (Long64_t i = 0; i < nEvt; ++i) {
      t->GetEntry(i);
      if (!te || !pe || te->GetEntries() == 0 || pe->GetEntries() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      auto *pidev = (AtPIDEvent *)pe->At(0);
      if (!ev || !pidev)
         continue;
      std::map<int, AtFittedTrack *> fmap;
      for (auto &ft : ev->GetFittedTracks())
         if (ft)
            fmap[ft->GetTrackID()] = ft.get();
      nFit += fmap.size();
      for (auto &sr : pidev->GetSpyral()) {
         if (!sr.valid || sr.brho <= 0)
            continue;
         auto it = fmap.find(sr.trackID);
         if (it == fmap.end())
            continue;
         auto *ft = it->second;
         const auto &m = ft->GetTrackMetadata();
         if (!m)
            continue;
         const double ndf = m->GetNdf(), chi2 = m->GetChi2();
         Rec r;
         r.keFit = ft->GetKinematics().kineticEnergy;
         r.keXtr = ft->GetKinematicsXtr().kineticEnergy;
         r.keBrho = keFromBrho(sr.brho);
         r.chi2n = ndf > 0 ? chi2 / ndf : 1e9;
         r.ndf = ndf;
         r.chi2 = chi2;
         r.matUsed = m->GetMatEffects();
         r.fellBack = m->GetMatEffectsFallback();
         if (r.keFit > 0 && r.keBrho > 0)
            out[i * 10000L + sr.trackID] = r;
      }
   }
   return out;
}

struct Stat {
   long n;
   double med, lo, hi;
};
Stat quant(std::vector<double> v)
{
   Stat s{(long)v.size(), 0, 0, 0};
   if (v.empty())
      return s;
   std::sort(v.begin(), v.end());
   auto at = [&](double q) { return v[std::min(v.size() - 1, (size_t)(q * v.size()))]; };
   s.med = at(0.50);
   s.lo = at(0.16);
   s.hi = at(0.84);
   return s;
}
void row(const char *lbl, Stat s)
{
   printf("  %-26s n=%-6ld  median=%7.4f   16-84%% = [%6.4f, %6.4f]  halfwidth=%6.4f\n", lbl, s.n, s.med, s.lo, s.hi,
          0.5 * (s.hi - s.lo));
}

} // namespace

void ab_catima_cmp(TString fOff, TString fOn)
{
   gSystem->Load("libAtReconstruction.so");
   long nEvtA = 0, nFitA = 0, nEvtB = 0, nFitB = 0;
   auto A = harvest(fOff, nEvtA, nFitA); // CATIMA off
   auto B = harvest(fOn, nEvtB, nFitB);  // CATIMA on

   printf("\n=== yields =========================================================\n");
   printf("  CATIMA off : %ld events, %ld fitted tracks, %zu with a Brho match\n", nEvtA, nFitA, A.size());
   printf("  CATIMA on  : %ld events, %ld fitted tracks, %zu with a Brho match\n", nEvtB, nFitB, B.size());

   long onlyA = 0, onlyB = 0, both = 0;
   for (auto &kv : A)
      (B.count(kv.first) ? both : onlyA)++;
   for (auto &kv : B)
      if (!A.count(kv.first))
         ++onlyB;
   printf("  paired: %ld   off-only: %ld   on-only: %ld\n", both, onlyA, onlyB);

   long matA = 0, fbA = 0, matB = 0, fbB = 0;
   for (auto &kv : A) {
      matA += kv.second.matUsed;
      fbA += kv.second.fellBack;
   }
   for (auto &kv : B) {
      matB += kv.second.matUsed;
      fbB += kv.second.fellBack;
   }
   printf("\n=== provenance (must be 100%% material, 0 fallback) =================\n");
   printf("  CATIMA off : matEffects used %ld/%zu, fell back %ld\n", matA, A.size(), fbA);
   printf("  CATIMA on  : matEffects used %ld/%zu, fell back %ld\n", matB, B.size(), fbB);

   // COLLAPSED-FIT CENSUS. This has to come before any width is quoted: a track with ndf < 0
   // carries kinematics into the ntuple like any other, so if one arm produces far more of them
   // than the other, a width comparison over all pairs is measuring the failure rate, not the
   // resolution.
   long badA = 0, badB = 0, badPairA = 0, badPairB = 0, bothGood = 0;
   for (auto &kv : A) {
      if (!kv.second.good())
         ++badA;
      auto ib = B.find(kv.first);
      if (ib == B.end())
         continue;
      const bool ga = kv.second.good(), gb = ib->second.good();
      if (!ga)
         ++badPairA;
      if (!gb)
         ++badPairB;
      if (ga && gb)
         ++bothGood;
   }
   for (auto &kv : B)
      if (!kv.second.good())
         ++badB;
   printf("\n=== collapsed fits (ndf <= 0 or chi2 <= 0) ==========================\n");
   printf("  CATIMA off : %ld/%zu (%.1f%%)      of the %ld paired: %ld\n", badA, A.size(), 100.0 * badA / A.size(),
          both, badPairA);
   printf("  CATIMA on  : %ld/%zu (%.1f%%)      of the %ld paired: %ld\n", badB, B.size(), 100.0 * badB / B.size(),
          both, badPairB);
   printf("  pairs good in BOTH arms: %ld\n", bothGood);

   const char *lbl[3] = {"all", "KE_Brho < 3.5 MeV", "KE_Brho > 3.5 MeV"};

   // pass 0 = every pair (contaminated by the collapsed fits above)
   // pass 1 = only pairs whose fit converged in BOTH arms -- the control
   for (int pass = 0; pass < 2; ++pass) {
      std::vector<double> rA[3], rB[3], rXA[3], rXB[3], cA, cB;
      for (auto &kv : A) {
         auto ib = B.find(kv.first);
         if (ib == B.end())
            continue;
         const Rec &a = kv.second, &b = ib->second;
         if (pass == 1 && !(a.good() && b.good()))
            continue;
         const int reg = (a.keBrho < kSplit) ? 1 : 2;
         for (int r : {0, reg}) {
            rA[r].push_back(a.keFit / a.keBrho);
            rB[r].push_back(b.keFit / b.keBrho);
            rXA[r].push_back(a.keXtr / a.keBrho);
            rXB[r].push_back(b.keXtr / b.keBrho);
         }
         cA.push_back(a.chi2n);
         cB.push_back(b.chi2n);
      }

      printf("\n########## %s ##########\n",
             pass == 0 ? "ALL PAIRS (contaminated -- shown only for contrast)"
                       : "CONTROL: pairs that CONVERGED IN BOTH ARMS");
      printf("\n=== KE_fit / KE_Brho  (raw fit) ====================================\n");
      for (int r = 0; r < 3; ++r) {
         printf("-- %s\n", lbl[r]);
         row("CATIMA off", quant(rA[r]));
         row("CATIMA on ", quant(rB[r]));
         Stat sa = quant(rA[r]), sb = quant(rB[r]);
         if (sa.n && sa.hi > sa.lo)
            printf("  %-26s median %+7.4f   halfwidth %+7.4f  (%+.1f%%)\n", "delta (on - off)", sb.med - sa.med,
                   0.5 * ((sb.hi - sb.lo) - (sa.hi - sa.lo)),
                   100.0 * (0.5 * (sb.hi - sb.lo) / (0.5 * (sa.hi - sa.lo)) - 1.0));
      }

      printf("\n=== KE_xtr / KE_Brho  (back-extrapolated, what the analysis uses) ==\n");
      for (int r = 0; r < 3; ++r) {
         printf("-- %s\n", lbl[r]);
         row("CATIMA off", quant(rXA[r]));
         row("CATIMA on ", quant(rXB[r]));
      }

      printf("\n=== chi2/ndf =======================================================\n");
      row("CATIMA off", quant(cA));
      row("CATIMA on ", quant(cB));
   }
   printf("\n");
}
