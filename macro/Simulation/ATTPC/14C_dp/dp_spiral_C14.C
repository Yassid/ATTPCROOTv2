/// @file dp_spiral_C14.C
/// @brief The local radius along a backward (d,p) proton track, window by window.
///
///   root -b -q 'dp_spiral_C14.C()'
///
/// The pattern circle radius falls with the number of hits in the track -- 0.991 of truth at 170
/// hits, 0.944 at 704, 0.888 at 1471 -- and the hits are 100 % pure proton, so nothing is being
/// merged in. The remaining explanation is that these protons are spiralling far enough to lose a
/// large fraction of their energy inside the chamber: a 5 MeV proton in D2 at 300 torr loses about
/// 1.2 keV/mm, so over a 2 m spiral it sheds half its energy and its radius shrinks as sqrt(E),
/// to about 0.7 of the initial value. One circle through all of that returns an average, and the
/// seed momentum -- which is supposed to be the momentum AT THE VERTEX -- inherits it.
///
/// This fits circles to successive windows of hits, ordered from the vertex, and reports each
/// window's radius against the truth radius at the vertex.

#include <algorithm>
#include <vector>

static bool sp_circle(const std::vector<double> &x, const std::vector<double> &y, double &R)
{
   const size_t n = x.size();
   if (n < 12) return false;
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
   double s = 0, s2 = 0;
   for (size_t i = 0; i < n; ++i) { double d = std::hypot(x[i] - cx, y[i] - cy); s += d; s2 += d * d; }
   R = s / N;
   double rms = std::sqrt(std::max(0.0, s2 / N - R * R));
   return R > 0 && rms < 0.35 * R; // reject windows that are not an arc at all
}
static double sp_q(std::vector<double> v, double p)
{
   if (v.size() < 5) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_spiral_C14(TString simFile = "/mnt/f/a1954_C14dp_hf/sims_b285/gs_s8001_sim.root",
                   TString recoFile = "/mnt/f/a1954_C14dp_matfx/ab_reco.root", Double_t bField = 2.85,
                   Int_t minHits = 600, Int_t nWin = 6)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   const double u = 931.49401, mp = 1.007825 * u, ZPAD = 1000.0;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   TFile *fr = TFile::Open(recoFile);
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   std::vector<double> win[16], geo, keAll;
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
      if (keT <= 0 || thT < 92 || thT > 140) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      const AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = &t; }
      if (!best || (int)nb < minHits) continue;
      std::vector<std::array<double, 3>> P;
      for (const auto &h : best->GetHitArray()) {
         auto q = h->GetPosition();
         P.push_back({q.X(), q.Y(), ZPAD - q.Z()});
      }
      std::sort(P.begin(), P.end(), [](const std::array<double, 3> &a, const std::array<double, 3> &c) {
         return a[2] > c[2];
      });
      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField);
      if (rTrue <= 0) continue;
      size_t w = P.size() / nWin;
      if (w < 12) continue;
      for (int j = 0; j < nWin; ++j) {
         std::vector<double> x, y;
         for (size_t k = j * w; k < (j + 1) * w && k < P.size(); ++k) { x.push_back(P[k][0]); y.push_back(P[k][1]); }
         double R;
         if (sp_circle(x, y, R)) win[j].push_back(R / rTrue);
      }
      geo.push_back(best->GetGeoRadius() / rTrue);
      keAll.push_back(keT);
   }

   printf("\n  Local circle radius / TRUE radius at the vertex, along backward tracks with >= %d hits\n", minHits);
   printf("  (%zu tracks, median KE %.1f MeV; windows are equal slices of the hit list from the vertex)\n\n",
          geo.size(), sp_q(keAll, .5));
   printf("  %-10s", "window");
   for (int j = 0; j < nWin; ++j) printf("  %4d/%d ", j + 1, nWin);
   printf("   | whole-track GeoRadius\n");
   printf("  %-10s", "R/R_true");
   for (int j = 0; j < nWin; ++j) {
      double v = sp_q(win[j], .5);
      if (std::isnan(v)) printf("  %6s ", "-"); else printf("  %6.3f ", v);
   }
   printf("   | %8.3f\n", sp_q(geo, .5));
   printf("  %-10s", "n");
   for (int j = 0; j < nWin; ++j) printf("  %6zu ", win[j].size());
   printf("\n\n  A falling sequence is the spiral tightening as the proton slows. The first window is the\n"
          "  one that carries the momentum the analysis wants.\n\n");
   fs->Close(); fr->Close();
}
