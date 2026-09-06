/// @file order_helix_dp.C
/// @brief INTRINSIC cluster ordering: sort by the UNWRAPPED azimuth of the helix.
///
/// WHY NEITHER EXISTING ORDERING WORKS (both measured on 800 long a1975 (d,p) tracks):
///   * AtGenfitter's z sort -- 90.6 % of tracks get a jump. At the jump the two clusters differ by
///     0.90 mm in z, below one time bucket (1.84 mm), so there is no z left to sort on and it picks
///     an arbitrary transverse neighbour; 15 % of the time on the opposite side of the circle.
///   * AtPRA's greedy 3D nearest-neighbour walk -- 88.9 % get a jump, median largest 124 mm. It
///     strands clusters and returns for them later (event 604: a 315 mm jump at index 218/244).
///
/// THE PARAMETER THAT ORDERS A HELIX IS ITS OWN PHASE. With axis along z,
///     z = z0 + (h/2pi) * Phi,     Phi = phi + 2*pi*k   (phi wrapped, k the turn number)
/// so for the CORRECT pitch h the quantity  u = z - (h/2pi)*phi  equals z0 + h*k: it collapses
/// onto a lattice of spacing h. That is a 1-D search needing NO prior ordering -- scan h, keep the
/// one whose u mod h is most concentrated, then k = round((u - z0)/h) and sort by Phi.
///
/// Nothing here uses a seed, a greedy choice, or the z resolution, which is exactly what the two
/// existing methods depend on.
///
/// KNOWN APPROXIMATION: a decelerating spiral has a SHRINKING radius, and its pitch shrinks with
/// it (the pitch angle is preserved by the magnetic force, so h/R stays roughly constant while
/// |p| falls). A single global h is therefore only an approximation; it is good enough as long as
/// the accumulated phase error stays below pi across the track, which is what the test below
/// measures rather than assumes.
#include <cmath>
#include <vector>
#include <numeric>

/// Kasa circle fit on the transverse projection.
static bool ohCircle(const std::vector<ROOT::Math::XYZPoint> &P, double &cx, double &cy, double &R)
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
   for (auto &p:P) R += std::hypot(p.X()-cx, p.Y()-cy);
   R/=m;
   return R>1.0;
}

/// Concentration of u = z - (h/2pi)*phi modulo h, as a resultant length in [0,1]; 1 = perfect.
static double ohConcentration(const std::vector<double> &z, const std::vector<double> &phi, double h)
{
   const double ah = std::fabs(h);
   if (ah < 1e-6) return 0;
   double sx=0, sy=0;
   const int n = z.size();
   for (int i=0;i<n;++i){
      double u = z[i] - (h/(2*M_PI))*phi[i];
      double a = 2*M_PI * (u/ah - std::floor(u/ah));
      sx += std::cos(a); sy += std::sin(a);
   }
   return std::hypot(sx,sy)/n;
}

/// The ordering. Returns indices sorted along the track, or an empty vector if it cannot decide.
/// `hOut` receives the fitted pitch (mm per turn, signed).
static std::vector<int> orderHelix(const std::vector<ROOT::Math::XYZPoint> &P, double &hOut,
                                   double *ROut=nullptr, double *concOut=nullptr)
{
   hOut = 0;
   const int n = P.size();
   std::vector<int> out;
   if (n < 8) return out;
   double cx,cy,R;
   if (!ohCircle(P,cx,cy,R)) return out;
   if (ROut) *ROut = R;
   std::vector<double> phi(n), z(n);
   double zlo=1e30, zhi=-1e30;
   for (int i=0;i<n;++i){
      phi[i] = std::atan2(P[i].Y()-cy, P[i].X()-cx);
      z[i]   = P[i].Z();
      zlo=std::min(zlo,z[i]); zhi=std::max(zhi,z[i]);
   }
   const double zspan = zhi-zlo;
   if (zspan < 1.0) return out;   // a track flat in z has no pitch to find
   // Scan the number of turns rather than h directly: h = zspan / turns, both signs.
   double bestC=-1, bestH=0;
   for (int s=0;s<2;++s){
      const double sign = s? -1.0 : 1.0;
      for (double turns=0.05; turns<=30.0; turns*=1.02){
         const double h = sign*zspan/turns;
         const double c = ohConcentration(z,phi,h);
         if (c>bestC){ bestC=c; bestH=h; }
      }
   }
   if (concOut) *concOut = bestC;
   if (bestC < 0.5) return out;    // no pitch describes these points: refuse rather than guess
   hOut = bestH;
   // turn numbers about the most populated lattice site
   const double ah=std::fabs(bestH);
   double sx=0, sy=0;
   for (int i=0;i<n;++i){ double u=z[i]-(bestH/(2*M_PI))*phi[i];
      double a=2*M_PI*(u/ah-std::floor(u/ah)); sx+=std::cos(a); sy+=std::sin(a); }
   const double a0 = std::atan2(sy,sx);
   const double u0 = ah*(a0<0? a0+2*M_PI : a0)/(2*M_PI);
   std::vector<double> Phi(n);
   for (int i=0;i<n;++i){
      double u=z[i]-(bestH/(2*M_PI))*phi[i];
      double k=std::round((u-u0)/bestH);
      Phi[i]=phi[i]+2*M_PI*k;
   }
   out.resize(n);
   std::iota(out.begin(),out.end(),0);
   std::sort(out.begin(),out.end(),[&](int a,int b){ return Phi[a]<Phi[b]; });
   return out;
}
