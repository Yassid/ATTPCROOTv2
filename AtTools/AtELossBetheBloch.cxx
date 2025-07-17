#include "AtELossBetheBloch.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace AtTools {
void AtELossBetheBloch::SetMaterial(int mat_A, int mat_Z, double density)
{
   fMat_A = mat_A;
   fMat_Z = mat_Z;
   SetDensity(density);
}

double AtELossBetheBloch::GetdEdx(double energy) const
{
   double gamma = energy / fPart_mass + 1; // This is KE not total energy
   double gammaSquare = gamma * gamma;
   double betaSquare = 1 - 1 / gammaSquare;
   std::cout << "Gamma: " << gamma << " Beta^2: " << betaSquare << " Energy: " << energy << " Part_mass: " << fPart_mass
             << std::endl;
   static const double betaGammaMin(0.05);
   if (betaSquare * gammaSquare < betaGammaMin * betaGammaMin) {
      throw std::invalid_argument("beta*gamma (" + std::to_string(std::sqrt(betaSquare * gammaSquare)) +
                                  ") < 0.05, Bethe-Bloch implementation not valid anymore!");
   }

   // calc dEdx_, also needed in noiseBetheBloch!
   double density = fDensity / 1000.; // Scale density from mg/cm^3 to g/cm^3
   double result(0.307075 * fMat_Z / fMat_A * density / betaSquare * fPart_q * fPart_q);
   double massRatio(fM_e / fMat_A);
   double argument(gammaSquare * betaSquare * fPart_mass * 1.E3 * 2. /
                   ((1.E-6 * fMat_A) * std::sqrt(1. + 2. * gamma * massRatio + massRatio * massRatio)));
   result *= std::log(argument) - betaSquare; // Bethe-Bloch [MeV/cm]
   if (result < 0.) {
      result = 0;
   }

   return result;
}

double AtELossBetheBloch::GetEnergyLoss(double energyIni, double distance) const
{
   // using fourth order Runge Kutta
   // k1 = f(t0,y0)
   // k2 = f(t0 + h/2, y0 + h/2 * k1)
   // k3 = f(t0 + h/2, y0 + h/2 * k2)
   // k4 = f(t0 + h,   y0 + h   * k3)

   // This means in our case:
   // dEdx1 = dEdx(x0,       E0)
   // dEdx2 = dEdx(x0 + h/2, E1); E1 = E0 + h/2 * dEdx1
   // dEdx3 = dEdx(x0 + h/2, E2); E2 = E0 + h/2 * dEdx2
   // dEdx4 = dEdx(x0 + h,   E3); E3 = E0 + h   * dEdx3
   return 0;
}

double AtELossBetheBloch::GetRange(double energyIni, double energyFin) const
{
   throw std::runtime_error("GetRange not implemented for Bethe-Bloch model");
   return 0;
}
double AtELossBetheBloch::GetEnergy(double energyIni, double distance) const
{
   return energyIni - GetEnergyLoss(energyIni, distance);
}

} // namespace AtTools
