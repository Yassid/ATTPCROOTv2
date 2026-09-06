/// @file theta_order_test.C
/// @brief Does ORDERING THE RAW HITS change the arc-length theta that AtPRA fits?
///
/// AtPRA::SetTrackInitialParameters builds the arc length from an UNWRAPPED azimuth, and the
/// unwrapping accumulates -Y-axis crossings IN HIT-ARRAY ORDER:
///     if (posOnCircle.X() < 0 && lastYSign != currYSign) numYCross -= currYSign;
///     angleHit += 2*M_PI*numYCross;
/// That is only meaningful if consecutive entries of the hit array are consecutive ALONG THE
/// TRACK. They are not: AtTrackFinderHDBSCAN fills the track from its cluster's point-index list
/// and calls SetTrackInitialParameters immediately, before any ordering. Measured on
/// run_0020 reco: the stored hit path is up to 2.8x the nearest-neighbour path, i.e. badly
/// scrambled. OrderClustersAlongTrack does NOT fix this -- it reorders the CLUSTER array, and the
/// arc-length fit reads the HIT array.
///
/// This reproduces the calculation offline both ways so the effect can be measured without a
/// re-reco. The line fit here is least squares rather than AtPRA's RANSAC, so the absolute values
/// differ a little from the stored GeoTheta; what is being compared is ordered vs unordered with
/// everything else identical.
#include <algorithm>
#include <cmath>
#include <vector>

static int tosign(double v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

/// circle fit (Kasa) on the transverse projection
static bool fitCirc(const std::vector<ROOT::Math::XYZPoint> &P, double &cx, double &cy, double &R)
{
   const int m = P.size();
   if (m < 5) return false;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
   for (auto &p : P) { double x=p.X(), y=p.Y();
      Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y; }
   double C=m*Sxx-Sx*Sx, D=m*Sxy-Sx*Sy, E=m*Sxxx+m*Sxyy-(Sxx+Syy)*Sx;
   double G=m*Syy-Sy*Sy, H=m*Sxxy+m*Syyy-(Sxx+Syy)*Sy;
   double den=2*(C*G-D*D);
   if (std::fabs(den)<1e-9) return false;
   cx=(E*G-D*H)/den; cy=(C*H-D*E)/den; R=0;
   for (auto &p : P) R += std::hypot(p.X()-cx, p.Y()-cy);
   R/=m;
   return R>0;
}

/// AtPRA's arc-length theta, verbatim in structure, on whatever order `P` is in
static double arcTheta(const std::vector<ROOT::Math::XYZPoint> &P)
{
   double cx, cy, R;
   if (!fitCirc(P, cx, cy, R)) return std::nan("");
   ROOT::Math::XYZPoint c(cx, cy, 0);
   double refAng = std::atan2(P[0].Y()-cy, P[0].X()-cx);
   int numYCross = 0, lastYSign = tosign(P[0].Y()-cy);
   std::vector<double> arc, z;
   for (size_t i = 0; i < P.size(); ++i) {
      double dx = P[i].X()-cx, dy = P[i].Y()-cy;
      double a = std::atan2(dy, dx);
      int currYSign = tosign(dy);
      if (dx < 0 && lastYSign != currYSign) numYCross -= currYSign;
      lastYSign = currYSign;
      a += 2*M_PI*numYCross;
      arc.push_back(R*(refAng - a));
      z.push_back(P[i].Z());
   }
   // least-squares line z = a + b*arc, then AtPRA's angle rule on the direction vector
   const int n = arc.size();
   double Sx=0,Sy=0,Sxx=0,Sxy=0;
   for (int i=0;i<n;++i){ Sx+=arc[i]; Sy+=z[i]; Sxx+=arc[i]*arc[i]; Sxy+=arc[i]*z[i]; }
   double den = n*Sxx - Sx*Sx;
   if (std::fabs(den) < 1e-9) return std::nan("");
   double b = (n*Sxy - Sx*Sy)/den;
   ROOT::Math::XYVector dir = ROOT::Math::XYVector(1.0, b).Unit();
   int sign = (dir.X()*dir.Y() < 0) ? -1 : 1;
   return std::acos(sign*std::fabs(dir.Y())) * TMath::RadToDeg();
}

void theta_order_test(TString reco = "/mnt/f/a1975/reco_d2_dv1104/run_0020_multifit_reco.root",
                      int maxTracks = 400, int minHits = 30)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto f = TFile::Open(reco);
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);
   int n = 0, flipped = 0, nanNow = 0, nanFix = 0, big = 0;
   double sumAbs = 0;
   printf("\n  %8s %6s %10s %10s %10s %9s\n","entry","nhits","GeoTheta","theta(raw)","theta(ord)","d(ord-raw)");
   for (Long64_t i = 0; i < t->GetEntries() && n < maxTracks; ++i) {
      t->GetEntry(i);
      if (!pe || pe->GetEntries()==0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      for (auto &tr : ev->GetTrackCand()) {
         auto &hits = tr.GetHitArray();
         if ((int)hits.size() < minHits) continue;
         std::vector<ROOT::Math::XYZPoint> S;
         for (auto &h : hits) S.push_back(h->GetPosition());
         // ORDERED copy: nearest-neighbour walk from the highest-z hit, the same convention
         // AtPRA::OrderClustersAlongTrack uses on clusters
         std::vector<ROOT::Math::XYZPoint> W;
         std::vector<char> used(S.size(), 0);
         int cur = 0;
         for (size_t k = 1; k < S.size(); ++k) if (S[k].Z() > S[cur].Z()) cur = k;
         used[cur]=1; W.push_back(S[cur]);
         for (size_t k = 1; k < S.size(); ++k) {
            int nx=-1; double bd=1e30;
            for (size_t j = 0; j < S.size(); ++j) { if (used[j]) continue;
               double d=(S[j]-W.back()).R(); if (d<bd){bd=d;nx=j;} }
            if (nx<0) break; used[nx]=1; W.push_back(S[nx]);
         }
         double tRaw = arcTheta(S), tOrd = arcTheta(W);
         double geo = tr.GetGeoTheta()*TMath::RadToDeg();
         if (std::isnan(geo)) ++nanNow;
         if (std::isnan(tOrd)) ++nanFix;
         if (!std::isnan(tRaw) && !std::isnan(tOrd)) {
            double d = tOrd - tRaw;
            sumAbs += std::fabs(d);
            if ((tRaw-90)*(tOrd-90) < 0) ++flipped;   // crosses the forward/backward boundary
            if (std::fabs(d) > 10) ++big;
            if (n < 12)
               printf("  %8lld %6zu %10.1f %10.1f %10.1f %9.1f%s\n", i, hits.size(), geo, tRaw, tOrd, d,
                      (tRaw-90)*(tOrd-90) < 0 ? "   <-- crosses 90 deg" : "");
         }
         ++n;
         break;
      }
   }
   printf("\n  %d tracks with >= %d hits\n", n, minHits);
   printf("  mean |theta(ordered) - theta(raw)| = %.1f deg\n", n ? sumAbs/n : 0.0);
   printf("  |change| > 10 deg            : %d  (%.0f%%)\n", big, 100.0*big/std::max(n,1));
   printf("  CROSSES the 90 deg boundary  : %d  (%.0f%%)  <- these flip forward/backward\n",
          flipped, 100.0*flipped/std::max(n,1));
   printf("  stored GeoTheta is NaN on    : %d\n", nanNow);
}
