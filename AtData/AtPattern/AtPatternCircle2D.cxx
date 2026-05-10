#include "AtPatternCircle2D.h"

#include "AtPadPlaneElement.h"

#include <Math/Point2D.h>
#include <Math/Point2Dfwd.h> // for XYPoint
#include <Math/Point3D.h>
#include <Math/Vector2D.h>    // for DisplacementVector2D
#include <Math/Vector3D.h>    // for DisplacementVector3D
#include <Math/Vector3Dfwd.h> // for RhoZPhiVector
#include <TEveLine.h>

#include <algorithm> // for max
#include <cmath>     // for fabs, isfinite, sqrt
#include <memory>    // for allocator_traits<>::value_type

class TEveElement;
using XYPoint = ROOT::Math::XYPoint;
using namespace AtPatterns;

AtPatternCircle2D::AtPatternCircle2D() : AtPattern(3) {}
TEveElement *AtPatternCircle2D::GetEveElement() const
{
   return AtPattern::GetEveLine(0, 2 * M_PI, 1000);
}

std::unique_ptr<AtPadPlaneElement> AtPatterns::AtPatternCircle2D::GetPadPlaneProjection() const
{
   return std::make_unique<AtPadPlaneCircle>(GetCenter().X(), GetCenter().Y(), GetRadius());
}

void AtPatternCircle2D::DefinePattern(const std::vector<XYZPoint> &points)
{
   XYPoint p1 = {points[0].X(), points[0].Y()};
   XYPoint p2 = {points[1].X(), points[1].Y()};
   XYPoint p3 = {points[2].X(), points[2].Y()};

   double a = 2 * (p2 - p1).X();
   double b = 2 * (p2 - p1).Y();
   double c = p2.Mag2() - p1.Mag2();
   double d = 2 * (p3 - p2).X();
   double e = 2 * (p3 - p2).Y();
   double f = p3.Mag2() - p2.Mag2();

   XYPoint center((b * f - e * c), (d * c - a * f));
   center /= (b * d - e * a);

   fPatternPar.clear();
   fPatternPar.push_back(center.X());
   fPatternPar.push_back(center.Y());
   fPatternPar.push_back((center - p1).R());
}

Double_t AtPatternCircle2D::DistanceToPattern(const XYZPoint &point) const
{
   auto pointToCenter = point - GetCenter();
   return std::abs(pointToCenter.Rho() - GetRadius());
}

XYZPoint AtPatternCircle2D::ClosestPointOnPattern(const XYZPoint &point) const
{
   auto pointToCenter = point - GetCenter();
   return GetPointAt(pointToCenter.Phi());
}

XYZPoint AtPatternCircle2D::GetPointAt(double theta) const
{
   return XYZPoint(GetCenter() + ROOT::Math::RhoZPhiVector(GetRadius(), 0, theta));
}

void AtPatternCircle2D::FitPattern(const std::vector<XYZPoint> &points, const std::vector<double> &charge)
{

   //  2D circle fit, due to Taubin, based on the journal article
   // G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   //             Space Curves Defined By Implicit Equations, With
   //             Applications To Edge And Range Image Segmentation",
   //             IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   //----- adapted from: https://people.cas.uab.edu/~mosya/cl/CPPcircle.html

   int iter, IterMAX = 99;
   int Niliers = 0;
   double Xi, Yi, Zi;
   double Xm = 0, Ym = 0, Zm = 0;
   double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
   double A0, A1, A2, A22, A3, A33;
   double Dy, xnew, x, ynew, y;
   double DET, Xcenter, Ycenter;
   double Q = 0;
   double a, b, r;

   bool doChargeWeight = points.size() == charge.size();

   for (int i = 0; i < points.size(); ++i) {
      const auto hitQ = doChargeWeight ? charge[i] : 1;
      const auto &pos = points[i];
      Q += hitQ / 10.;
      Xm += pos.X() * hitQ / 10.;
      Ym += pos.Y() * hitQ / 10.;
      // Zm += pos.Z() * hitQ / 10.;
   }

   Xm /= Q;
   Ym /= Q;
   // Zm /= Q;

   Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;

   for (const auto &pos : points) {
      Xi = pos.X() - Xm; //  centered x-coordinates
      Yi = pos.Y() - Ym; //  centered y-coordinates
      Zi = Xi * Xi + Yi * Yi;

      Mxy += Xi * Yi;
      Mxx += Xi * Xi;
      Myy += Yi * Yi;
      Mxz += Xi * Zi;
      Myz += Yi * Zi;
      Mzz += Zi * Zi;
   }

   Mxx /= points.size();
   Myy /= points.size();
   Mxy /= points.size();
   Mxz /= points.size();
   Myz /= points.size();
   Mzz /= points.size();

   //      computing coefficients of the characteristic polynomial

   Mz = Mxx + Myy;
   Cov_xy = Mxx * Myy - Mxy * Mxy;
   Var_z = Mzz - Mz * Mz;
   A3 = 4.0 * Mz;
   A2 = -3.0 * Mz * Mz - Mzz;
   A1 = Var_z * Mz + 4.0 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
   A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
   A22 = A2 + A2;
   A33 = A3 + A3 + A3;

   //    finding the root of the characteristic polynomial
   //    using Newton's method starting at x=0
   //     (it is guaranteed to converge to the right root)

   for (x = 0., y = A0, iter = 0; iter < IterMAX; iter++) // usually, 4-6 iterations are enough
   {
      Dy = A1 + x * (A22 + A33 * x);
      xnew = x - y / Dy;
      if ((xnew == x) || (!std::isfinite(xnew)))
         break;
      ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (std::abs(ynew) >= std::abs(y))
         break;
      x = xnew;
      y = ynew;
   }

   //       computing paramters of the fitting circle

   DET = x * x - x * Mz + Cov_xy;
   Xcenter = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2.0;
   Ycenter = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2.0;

   //       assembling the output

   a = Xcenter + Xm;
   b = Ycenter + Ym;
   r = sqrt(Xcenter * Xcenter + Ycenter * Ycenter + Mz);

   // ----------------------------------------------------------------------
   // Geometric refinement on top of Taubin's algebraic fit.
   //
   // Taubin minimises algebraic residuals (xi² + yi² − 2A·xi − 2B·yi − C)²
   // weighted by curvature, which gives a closed-form solution but biases
   // toward smaller R for weakly-curved arcs (large R, small sagitta).
   // For pi-TPC MIPs at B=0.5 T (R ≈ 400 mm, sagitta ~3 mm over 100 mm
   // chord) Taubin gives σ(R)/R ≈ 10% even though the geometric floor is
   // ~3%. A few Gauss-Newton iterations of the orthogonal-distance fit
   //   minimise Σ wi · (√[(xi−xc)² + (yi−yc)²] − R)²
   // close most of that gap. For tightly-curved arcs (16C+pp protons at
   // 2.85 T) the GN refit converges in 1-2 iterations to nearly the same
   // values Taubin returns, so it's safe across physics.
   const int kGNMaxIter = 10;
   const double kGNTol = 1e-6; // relative parameter change to declare convergence
   for (int gnIter = 0; gnIter < kGNMaxIter; ++gnIter) {
      // Build normal equations for [Δxc, Δyc, ΔR]: J^T W J · Δp = -J^T W r.
      double H[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      double g[3] = {0, 0, 0};
      double sumW = 0;
      bool degenerate = false;
      for (size_t i = 0; i < points.size(); ++i) {
         const double w = doChargeWeight ? std::max(0.0, charge[i] / 10.0) : 1.0;
         const double dx = points[i].X() - a;
         const double dy = points[i].Y() - b;
         const double d = std::sqrt(dx * dx + dy * dy);
         if (d < 1e-9) { degenerate = true; break; }
         const double res = d - r;
         // Jacobian row: ∂res/∂xc = −dx/d, ∂res/∂yc = −dy/d, ∂res/∂R = −1
         const double jx = -dx / d;
         const double jy = -dy / d;
         const double jR = -1.0;
         H[0][0] += w * jx * jx;
         H[0][1] += w * jx * jy;
         H[0][2] += w * jx * jR;
         H[1][1] += w * jy * jy;
         H[1][2] += w * jy * jR;
         H[2][2] += w * jR * jR;
         g[0] -= w * jx * res;
         g[1] -= w * jy * res;
         g[2] -= w * jR * res;
         sumW += w;
      }
      if (degenerate || sumW <= 0) break;
      H[1][0] = H[0][1];
      H[2][0] = H[0][2];
      H[2][1] = H[1][2];

      // Solve 3x3 by Cramer (symmetric positive-definite expected).
      const double D = H[0][0] * (H[1][1] * H[2][2] - H[1][2] * H[2][1])
                       - H[0][1] * (H[1][0] * H[2][2] - H[1][2] * H[2][0])
                       + H[0][2] * (H[1][0] * H[2][1] - H[1][1] * H[2][0]);
      if (std::abs(D) < 1e-30) break;
      const double dxc = (g[0] * (H[1][1] * H[2][2] - H[1][2] * H[2][1])
                          - H[0][1] * (g[1] * H[2][2] - H[1][2] * g[2])
                          + H[0][2] * (g[1] * H[2][1] - H[1][1] * g[2])) / D;
      const double dyc = (H[0][0] * (g[1] * H[2][2] - H[1][2] * g[2])
                          - g[0] * (H[1][0] * H[2][2] - H[1][2] * H[2][0])
                          + H[0][2] * (H[1][0] * g[2] - g[1] * H[2][0])) / D;
      const double dR = (H[0][0] * (H[1][1] * g[2] - g[1] * H[2][1])
                         - H[0][1] * (H[1][0] * g[2] - g[1] * H[2][0])
                         + g[0] * (H[1][0] * H[2][1] - H[1][1] * H[2][0])) / D;

      a += dxc;
      b += dyc;
      r += dR;
      const double rel = (std::abs(dxc) + std::abs(dyc) + std::abs(dR))
                         / std::max(1.0, std::abs(r));
      if (rel < kGNTol) break;
      if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(r)) || r <= 0) {
         // GN diverged — restore Taubin solution and bail.
         a = Xcenter + Xm;
         b = Ycenter + Ym;
         r = std::sqrt(Xcenter * Xcenter + Ycenter * Ycenter + Mz);
         break;
      }
   }

   fPatternPar.clear();
   fPatternPar.push_back(a);
   fPatternPar.push_back(b);
   fPatternPar.push_back(r);
   fChi2 = fabs(Mz - r * r); // TODO: Is this a chi2 value?

   //  return fabs(Mz - r * r);
}
