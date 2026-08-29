/// @file dp_arc_diag_C14.C
/// @brief Does the pattern circle under-estimate a backward track's radius because the track is a
/// SHRINKING spiral?
///
///   root -b -q 'dp_arc_diag_C14.C()'
///
/// A backward (d,p) proton is slow (3-6 MeV) with a metres-long range in D2, so it curls through
/// several turns inside the chamber, losing energy as it goes: the projected hits lie on a spiral
/// of DECREASING radius, not on a circle. One circle fitted to all of it returns something between
/// the first and last radius, which is systematically below the radius at the vertex -- the only
/// one that carries the initial momentum.
///
/// The test: fit the circle to the first f of the hits, counted from the VERTEX end, and compare
/// with the truth transverse radius. If the first-arc radius is right while the whole-track radius
/// is 5-8 % low, the diagnosis is settled and the fix is to seed from the first arc.

#include <algorithm>
#include <vector>

/// algebraic (Kasa) circle fit; returns false if degenerate
static bool arc_circle(const std::vector<double> &x, const std::vector<double> &y, double &R)
{
   const size_t n = x.size();
   if (n < 5) return false;
   double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
   for (size_t i = 0; i < n; ++i) {
      double z = x[i] * x[i] + y[i] * y[i];
      Sx += x[i]; Sy += y[i]; Sxx += x[i] * x[i]; Syy += y[i] * y[i];
      Sxy += x[i] * y[i]; Sxz += x[i] * z; Syz += y[i] * z; Sz += z;
   }
   double N = (double)n;
   double a11 = 2 * (Sxx - Sx * Sx / N), a12 = 2 * (Sxy - Sx * Sy / N), a22 = 2 * (Syy - Sy * Sy / N);
   double b1 = Sxz - Sx * Sz / N, b2 = Syz - Sy * Sz / N;
   double det = a11 * a22 - a12 * a12;
   if (std::abs(det) < 1e-9) return false;
   double cx = (b1 * a22 - b2 * a12) / det, cy = (a11 * b2 - a12 * b1) / det;
   double s = 0;
   for (size_t i = 0; i < n; ++i) s += std::hypot(x[i] - cx, y[i] - cy);
   R = s / N;
   return R > 0;
}
static double arc_q(std::vector<double> v, double p)
{
   if (v.size() < 8) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_arc_diag_C14(TString simFile = "/mnt/f/a1954_C14dp_hf/sims_b285/gs_s8001_sim.root",
                     TString recoFile = "/mnt/f/a1954_C14dp_matfx/ab_reco.root", Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   const double u = 931.49401, mp = 1.007825 * u;
   const double ZPAD = 1000.0;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   TFile *fr = TFile::Open(recoFile);
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   const int NF = 5;
   const double frac[NF] = {0.10, 0.20, 0.35, 0.60, 1.00};
   const int NB = 3;
   const double lo[NB] = {60, 95, 115}, hi[NB] = {85, 115, 150};
   std::vector<double> rr[NB][NF];
   std::vector<double> turns[NB], hits[NB];

   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pm = std::sqrt(px * px + py * py + pz * pz);
         if (pm <= 0) continue;
         keT = std::sqrt(pm * pm + mp * mp) - mp;
         thT = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (keT <= 0) continue;
      int b = -1;
      for (int j = 0; j < NB; ++j) if (thT >= lo[j] && thT < hi[j]) b = j;
      if (b < 0) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      const AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = &t; }
      if (!best || nb < 30) continue;

      // lab-frame hits, ordered by z_lab; the vertex end is the HIGHEST z_lab for a backward track
      std::vector<std::array<double, 3>> P;
      for (const auto &h : best->GetHitArray()) {
         auto q = h->GetPosition();
         P.push_back({q.X(), q.Y(), ZPAD - q.Z()});
      }
      const bool backward = (thT > 90);
      std::sort(P.begin(), P.end(), [&](const std::array<double, 3> &a, const std::array<double, 3> &c) {
         return backward ? (a[2] > c[2]) : (a[2] < c[2]);
      });

      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField);

      for (int f = 0; f < NF; ++f) {
         size_t m = std::max<size_t>(8, (size_t)(frac[f] * P.size()));
         m = std::min(m, P.size());
         std::vector<double> x, y;
         for (size_t k = 0; k < m; ++k) { x.push_back(P[k][0]); y.push_back(P[k][1]); }
         double R;
         if (arc_circle(x, y, R) && rTrue > 0) rr[b][f].push_back(R / rTrue);
      }
      // how many turns: path length in the pad plane divided by the circumference of the truth circle
      double path = 0;
      for (size_t k = 1; k < P.size(); ++k)
         path += std::hypot(P[k][0] - P[k - 1][0], P[k][1] - P[k - 1][1]);
      if (rTrue > 0) turns[b].push_back(path / (2 * TMath::Pi() * rTrue));
      hits[b].push_back((double)P.size());
   }

   printf("\n  Pattern circle radius / TRUE transverse radius, as a function of how much of the\n"
          "  track is used, counted from the vertex end. 14C(d,p) at %.2f T.\n\n", bField);
   printf("  %-12s %6s %7s %7s |", "theta_lab", "n", "<hits>", "<turns>");
   for (int f = 0; f < NF; ++f) printf("  first %3.0f%%", 100 * frac[f]);
   printf("\n");
   for (int b = 0; b < NB; ++b) {
      if (rr[b][0].size() < 8) { printf("  %3.0f-%-8.0f (too few)\n", lo[b], hi[b]); continue; }
      printf("  %3.0f-%-8.0f %6zu %7.0f %7.2f |", lo[b], hi[b], rr[b][0].size(), arc_q(hits[b], .5),
             arc_q(turns[b], .5));
      for (int f = 0; f < NF; ++f) printf("  %10.3f", arc_q(rr[b][f], .5));
      printf("\n");
   }
   printf("\n  If the first-arc value is ~1.00 and the whole-track value is well below it, the pattern\n"
          "  circle is averaging over a spiral that tightens as the proton slows, and the seed momentum\n"
          "  inherits that. <turns> is the projected path length over the truth circumference.\n\n");
   fs->Close(); fr->Close();
}
