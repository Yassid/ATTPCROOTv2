#ifndef ATELOSSBETHEBLOCH_H
#define ATELOSSBETHEBLOCH_H

#include "AtELossModel.h"

namespace AtTools {

/**
 * Class representing the energy loss of a particle through some material.
 * Derived classes can represent models from different sources (SRIM, etc).
 * Internal units are Mev/mm
 *
 * Based on a combination of Nabin Rijal's AtELossManager and the EnergyLoss class
 * (https://github.com/joshhooker/EnergyLossClass) which is released unter the MIT Licsense (copyright Joshua Hooker).
 */

class AtELossBetheBloch : public AtELossModel {
protected:
   static constexpr double fM_e = 0.51099895069e-3; // Mass of electron in eV/c^2

   double fPart_q;    // charge of the particle in e
   double fPart_mass; // Mass of particle in MeV/c^2
   int fMat_Z;        // Atomic number of the material
   int fMat_A;        // Mass number of the material

public:
   AtELossBetheBloch(double part_q, double part_mass, int mat_Z, int mat_A, double density)
      : AtELossModel(density), fPart_q(part_q), fPart_mass(part_mass), fMat_Z(mat_Z), fMat_A(mat_A)
   {
   }
   virtual ~AtELossBetheBloch() = default;

   void SetMaterial(int mat_A, int mat_Z, double density);

   /**
    * Get the stopping power in MeV/mm
    * @param energy Energy of the particle in MeV
    * @return Stopping power in MeV/mm
    */
   virtual double GetdEdx(double energy) const override;

   /**
    * Get the range of the particle in the material.
    */
   virtual double GetRange(double energyIni, double energyFin = 0) const override;

   /**
    * Get the energy loss over some distance (in mm).
    */
   virtual double GetEnergyLoss(double energyIni, double distance) const override;

   /**
    * Get the energy of particle after traveling some distance (in mm).
    * If the distance is negative, then returns the energy the particle had to
    * reach energyIni after distance.
    */
   virtual double GetEnergy(double energyIni, double distance) const override;

protected:
   double integrateRK4(double stepSize, double length);
};

} // namespace AtTools

#endif
