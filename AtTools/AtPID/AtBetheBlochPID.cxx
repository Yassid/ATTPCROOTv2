#include "AtBetheBlochPID.h"

#include "AtHit.h"
#include "AtTrack.h"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace AtTools {

void AtBetheBlochPID::AddSpecies(const std::string &name, double massMeV, int absPDG)
{
   fSpecies.push_back({name, massMeV, absPDG});
}

double AtBetheBlochPID::BetheBlochShape(double p_MeV, double massMeV, double I_MeV)
{
   if (!(p_MeV > 0) || !(massMeV > 0))
      return 1e-6;
   const double me = 0.510999;            // electron mass, MeV
   double bg = p_MeV / massMeV;           // beta*gamma
   double bg2 = bg * bg;
   double b2 = bg2 / (1.0 + bg2);         // beta^2
   double g2 = 1.0 + bg2;                 // gamma^2
   double Tmax = 2 * me * bg2 / (1 + 2 * std::sqrt(g2) * me / massMeV + (me / massMeV) * (me / massMeV));
   double lnArg = 2 * me * b2 * g2 * Tmax / (I_MeV * I_MeV);
   if (!(lnArg > 0))
      return 1e-6;
   double val = (1.0 / b2) * (0.5 * std::log(lnArg) - b2);
   return val > 0 ? val : 1e-6;
}

void AtBetheBlochPID::SetReferenceCalibration(double refRigidityMeV, double refDeDx, double refMassMeV)
{
   double shape = BetheBlochShape(refRigidityMeV, refMassMeV);
   if (shape > 0 && refDeDx > 0)
      fK = refDeDx / shape;
}

double AtBetheBlochPID::ExpectedDeDx(double rigidityMeV, double massMeV) const
{
   return fK * BetheBlochShape(rigidityMeV, massMeV);
}

int AtBetheBlochPID::Classify(double rigidityMeV, double dEdx) const
{
   if (fSpecies.empty() || !(rigidityMeV > 0) || !(dEdx > 0))
      return -1;
   int best = -1;
   double bestDist = 1e300;
   double lnMeas = std::log(dEdx);
   for (size_t i = 0; i < fSpecies.size(); ++i) {
      double exp = ExpectedDeDx(rigidityMeV, fSpecies[i].massMeV);
      double d = std::abs(lnMeas - std::log(exp)); // nearest in log space
      if (d < bestDist) {
         bestDist = d;
         best = static_cast<int>(i);
      }
   }
   return best;
}

bool AtBetheBlochPID::Observables(const AtTrack &track, double bField_T, double &rigidityMeV, double &dEdx,
                                  double trimHighFrac)
{
   double R = const_cast<AtTrack &>(track).GetGeoRadius(); // mm
   if (!(R > 0 && R < 1e5))
      return false;
   // transverse rigidity p_T = 0.2998 * B * R (mass-independent, |q|=1)
   rigidityMeV = 0.299792458 * bField_T * R;

   auto cen = const_cast<AtTrack &>(track).GetGeoCenter();
   const auto &hits = track.GetHitArray();
   if (hits.size() < 4)
      return false;
   // order hits azimuthally about the fitted circle centre, then per-hit dQ/dl
   std::vector<std::tuple<double, double, double, double, double>> s; // phi,x,y,z,q
   s.reserve(hits.size());
   for (const auto &h : hits) {
      const auto &p = h->GetPosition();
      s.emplace_back(std::atan2(p.Y() - cen.second, p.X() - cen.first), p.X(), p.Y(), p.Z(), h->GetCharge());
   }
   std::sort(s.begin(), s.end());
   std::vector<double> samples;
   samples.reserve(s.size());
   for (size_t i = 1; i < s.size(); ++i) {
      double dx = std::get<1>(s[i]) - std::get<1>(s[i - 1]);
      double dy = std::get<2>(s[i]) - std::get<2>(s[i - 1]);
      double dz = std::get<3>(s[i]) - std::get<3>(s[i - 1]);
      double dl = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (dl > 0.05)
         samples.push_back(std::get<4>(s[i]) / dl);
   }
   if (samples.size() < 3)
      return false;
   std::sort(samples.begin(), samples.end());
   size_t keep = static_cast<size_t>(std::ceil((1.0 - trimHighFrac) * samples.size())); // drop Landau tail
   if (keep < 1)
      keep = 1;
   double sum = 0;
   for (size_t i = 0; i < keep; ++i)
      sum += samples[i];
   dEdx = sum / keep;
   return dEdx > 0;
}

} // namespace AtTools
