#include "AtDirDeDxCleaner.h"

#include "AtHit.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace AtTools {
namespace DataCleaning {

HitCloud AtDirDeDxCleaner::CleanData(const HitCloud &hits)
{
   const int n = static_cast<int>(hits.size());
   HitCloud ret;
   if (n < fK + 1) { // too few hits to judge structure -> keep all
      for (const auto &h : hits)
         ret.push_back(std::make_unique<AtHit>(*h));
      return ret;
   }

   std::vector<double> X(n), Y(n), Z(n), Q(n);
   for (int i = 0; i < n; ++i) {
      auto p = hits[i]->GetPosition();
      X[i] = p.X();
      Y[i] = p.Y();
      Z[i] = p.Z();
      Q[i] = hits[i]->GetCharge();
   }

   // kNN indices (brute force; n is a few hundred per event)
   std::vector<std::vector<int>> nb(n);
   std::vector<std::pair<double, int>> dd(n);
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
         double dx = X[j] - X[i], dy = Y[j] - Y[i], dz = Z[j] - Z[i];
         dd[j] = {dx * dx + dy * dy + dz * dz, j};
      }
      std::nth_element(dd.begin(), dd.begin() + fK + 1, dd.end());
      std::sort(dd.begin(), dd.begin() + fK + 1);
      nb[i].reserve(fK);
      for (int m = 1; m <= fK; ++m) // skip self at m==0
         nb[i].push_back(dd[m].second);
   }

   // per-hit local tangent (principal axis via power iteration) + smoothed charge
   std::vector<double> tx(n), ty(n), tz(n), Qeff(n);
   std::vector<double> qbuf;
   qbuf.reserve(fK);
   for (int i = 0; i < n; ++i) {
      double cx = 0, cy = 0, cz = 0;
      for (int j : nb[i]) {
         cx += X[j];
         cy += Y[j];
         cz += Z[j];
      }
      cx /= fK;
      cy /= fK;
      cz /= fK;
      double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
      for (int j : nb[i]) {
         double ax = X[j] - cx, ay = Y[j] - cy, az = Z[j] - cz;
         c00 += ax * ax;
         c01 += ax * ay;
         c02 += ax * az;
         c11 += ay * ay;
         c12 += ay * az;
         c22 += az * az;
      }
      double vx = 1, vy = 1, vz = 1;
      for (int it = 0; it < 24; ++it) {
         double wx = c00 * vx + c01 * vy + c02 * vz;
         double wy = c01 * vx + c11 * vy + c12 * vz;
         double wz = c02 * vx + c12 * vy + c22 * vz;
         double nrm = std::sqrt(wx * wx + wy * wy + wz * wz);
         if (nrm < 1e-12)
            break;
         vx = wx / nrm;
         vy = wy / nrm;
         vz = wz / nrm;
      }
      tx[i] = vx;
      ty[i] = vy;
      tz[i] = vz;
      if (fSmoothQ) {
         qbuf.clear();
         for (int j : nb[i])
            qbuf.push_back(Q[j]);
         std::nth_element(qbuf.begin(), qbuf.begin() + qbuf.size() / 2, qbuf.end());
         Qeff[i] = qbuf[qbuf.size() / 2];
      } else {
         Qeff[i] = Q[i];
      }
   }

   // degree = number of direction- AND dE/dx-continuous incident edges
   std::vector<int> deg(n, 0);
   for (int i = 0; i < n; ++i) {
      for (int j : nb[i]) {
         double sx = X[j] - X[i], sy = Y[j] - Y[i], sz = Z[j] - Z[i];
         double dl = std::sqrt(sx * sx + sy * sy + sz * sz);
         if (dl > fRMax || dl < 1e-9)
            continue;
         double ux = sx / dl, uy = sy / dl, uz = sz / dl;
         double ai = std::fabs(ux * tx[i] + uy * ty[i] + uz * tz[i]);
         double aj = std::fabs(ux * tx[j] + uy * ty[j] + uz * tz[j]);
         double at = std::fabs(tx[i] * tx[j] + ty[i] * ty[j] + tz[i] * tz[j]);
         double qr = std::min(Qeff[i], Qeff[j]) / (std::max(Qeff[i], Qeff[j]) + 1e-9);
         if (ai > fCosSeg && aj > fCosSeg && at > fCosTan && qr > fQRatio) {
            deg[i]++;
            deg[j]++;
         }
      }
   }

   // keep a hit if it has enough continuous neighbours OR its charge exceeds the override
   // (high-dE/dx signal is kept even when isolated; noise is baseline-low so this costs ~0 noise)
   for (int i = 0; i < n; ++i)
      if (deg[i] >= fMinDeg || Q[i] > fQKeep)
         ret.push_back(std::make_unique<AtHit>(*hits[i]));
   return ret;
}

} // namespace DataCleaning
} // namespace AtTools
