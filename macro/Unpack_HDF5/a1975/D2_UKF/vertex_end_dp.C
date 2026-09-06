/// @file vertex_end_dp.C
/// @brief Which end of a track is the VERTEX end -- by EXTRAPOLATION to the beam axis, with no
///        fit and no charge/radius heuristic.
///
/// THE ARGUMENT (Yassid, 2026-09-05). Charge and curvature only CORRELATE with the direction.
/// The vertex is *defined* by the beam axis, so the direct test is: continue the trajectory
/// outward past each end and see which continuation reaches x = y = 0. The track is a helix, so
/// its transverse projection is a circle of radius R about a guiding centre C -- the circle's
/// closest approach to the beam axis is a single point P*, on the line from the origin through C:
///
///     P* = C * (|C| - R) / |C|          distance of closest approach  d* = | |C| - R |
///
/// The vertex end is then the one from which P* is reached FIRST when continuing away from the
/// track body. That is a pure geometric question about arc angles -- no fit, no energy loss, no
/// z convention, and it works on tracks whose fit failed, which is exactly where the existing
/// direction determination cannot be trusted.
///
/// WHY THE EXISTING DIRECTION CANNOT BE USED HERE. AtPRA::SetTrackInitialParameters gets theta
/// from a line fit in the (arc length, z) plane, and it builds the arc length by unwrapping the
/// azimuth IN HIT-ARRAY ORDER. AtTrackFinderHDBSCAN fills the hit array from its cluster's
/// point-index list and calls that immediately, before any ordering, so the unwrapping runs along
/// a scrambled path (measured on run_0020: stored hit path up to 2.8x the nearest-neighbour path).
/// Reordering the hits first moves theta by more than 10 deg on 40% of tracks and across the
/// 90 deg forward/backward boundary on 32%. OrderClustersAlongTrack does not help: it reorders the
/// CLUSTER array while the arc-length fit reads the HIT array.
///
///   root -b -q 'vertex_end_dp.C("list.txt","out.csv")'
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

struct VeInfo {
   bool ok{false};
   bool reverse{false};   ///< true = the LAST cluster is the vertex end, so reverse the array
   double dStar{0};       ///< closest approach of the fitted circle to the beam axis (mm)
   double arcFront{0};    ///< arc angle from the FIRST cluster outward to P* (rad)
   double arcBack{0};     ///< arc angle from the LAST cluster outward to P* (rad)
   double R{0}, cx{0}, cy{0};
   double spanRad{0};     ///< angular span covered by the track itself
};

static bool veFitCircle(const std::vector<ROOT::Math::XYZPoint> &P, double &cx, double &cy, double &R)
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
   return R > 1.0;
}

/// signed azimuth difference b - a, wrapped to (-pi, pi]
static double veWrap(double d) { while (d > M_PI) d -= 2*M_PI; while (d <= -M_PI) d += 2*M_PI; return d; }

VeInfo veVertexEnd(const std::vector<ROOT::Math::XYZPoint> &P)
{
   VeInfo I;
   if (!veFitCircle(P, I.cx, I.cy, I.R)) return I;
   const double cR = std::hypot(I.cx, I.cy);
   if (cR < 1e-6) return I;                       // guiding centre on the axis: P* undefined
   I.dStar = std::fabs(cR - I.R);
   // P*, the point of the circle closest to the beam axis
   const double px = I.cx * (cR - I.R) / cR, py = I.cy * (cR - I.R) / cR;
   const double aStar = std::atan2(py - I.cy, px - I.cx);

   // unwrap the track's own azimuth so its angular span and sense of travel are known
   std::vector<double> a(P.size());
   for (size_t i = 0; i < P.size(); ++i) a[i] = std::atan2(P[i].Y()-I.cy, P[i].X()-I.cx);
   double acc = a[0], span = 0;
   std::vector<double> u(P.size()); u[0] = acc;
   for (size_t i = 1; i < P.size(); ++i) { acc += veWrap(a[i]-a[i-1]); u[i] = acc; }
   span = u.back() - u.front();
   I.spanRad = std::fabs(span);
   const int sense = (span >= 0) ? +1 : -1;        // +1 = azimuth increases along the array

   // Continue OUTWARD past each end -- past the front means stepping backwards in azimuth,
   // past the back means stepping forwards -- and measure the angle to P* in that direction.
   auto outward = [&](double from, int dir) {
      double d = veWrap(aStar - from) * dir;
      if (d < 0) d += 2*M_PI;                      // always the outward-going angle
      return d;
   };
   I.arcFront = outward(u.front(), -sense);
   I.arcBack  = outward(u.back(),  +sense);
   // the vertex end is the one that reaches the beam axis sooner going outward
   I.reverse = (I.arcBack < I.arcFront);
   I.ok = true;
   return I;
}

void vertex_end_dp(TString listFile = "/mnt/f/a1975/caches/.explorer/refit_list.txt",
                   TString csvOut = "/tmp/vend.csv",
                   TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while (in >> r >> e >> ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit", e, ti});
   std::ofstream o(csvOut.Data());
   o << "run,entry,tid,ncl,ok,reverse,dstar,arcfront,arcback,span,R,ccen\n";
   TString open; TFile *f=nullptr; TTree *t=nullptr; TClonesArray *te=nullptr; long n=0;
   for (auto &j : jobs) {
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if (rt!=open) { if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); te=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&te); open=rt; }
      t->GetEntry(en);
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      auto trks=ev->GetTrackArray(); AtTrack *src=nullptr;
      for (auto &tr:trks) if (tr.GetTrackID()==id) src=&tr;
      if(!src) continue;
      auto *hc=src->GetHitClusterArray(); if(!hc||hc->size()<8) continue;
      std::vector<ROOT::Math::XYZPoint> P;
      for (auto &c : *hc) P.push_back(c.GetPosition());
      VeInfo I = veVertexEnd(P);
      o << rt << "," << en << "," << id << "," << P.size() << "," << (I.ok?1:0) << ","
        << (I.reverse?1:0) << "," << I.dStar << "," << I.arcFront << "," << I.arcBack << ","
        << I.spanRad << "," << I.R << "," << std::hypot(I.cx,I.cy) << "\n";
      ++n;
   }
   if (f) f->Close();
   o.close();
   printf("  wrote %s  (%ld tracks)\n", csvOut.Data(), n);
}
