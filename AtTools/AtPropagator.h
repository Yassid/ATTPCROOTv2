#ifndef ATPROPAGATOR_H
#define ATPROPAGATOR_H

#include "AtELossModel.h"

#include "Math/Vector3D.h"

namespace AtTools {

/**
 * @brief Class for propagating particles through a medium.
 *
 * This class is responsible for simulating the propagation of particles
 * through a medium, taking into account energy loss and other effects.
 * Uses an AtELossModel to calculate the energy loss and propagates particles
 * in the presence of electric and magnetic fields.
 *
 * Class is designed to be used with a single particle type. Create a new instance of the
 * class if the material or particle type changes.
 */
class AtPropagator {
protected:
   using XYZVector = ROOT::Math::XYZVector;
   XYZVector fEField{0, 0, 0}; // Electric field vector
   XYZVector fBField{0, 0, 0}; // Magnetic field vector

   const double fQ;                                 // Charge of the particle in Coulombs
   const double fMass;                              // Mass of the particle in MeV/c^2
   const std::unique_ptr<AtELossModel> fELossModel; // Energy loss model

   // Internal state variables for the propagator
   double fH = 1e-10;           /// Step size for propagation in s
   double fETol = 1e-4;         /// Energy tolerance for convergence when fixing energy loss
   double fStopTol = 0.01;      /// Maximum kinetic energy to consider the particle stopped
   double fScalingFactor = 1.0; /// Scaling factor for energy loss

   XYZVector fPos; // Current position of the particle in mm
   XYZVector fMom; // Current momentum of the particle in MeV/c

   static constexpr double fReltoSImom = 1.60218e-13 / 299792458; // Conversion factor from MeV/c to kg m/s (SI units)

public:
   AtPropagator(double charge, double mass, std::unique_ptr<AtELossModel> elossModel)
      : fQ(charge), fMass(mass), fELossModel(std::move(elossModel))
   {
   }
   /**
    * @brief Set the electric field (V/m)
    */
   void SetEField(const XYZVector &eField) { fEField = eField; }
   void SetEField(double ex, double ey, double ez) { fEField.SetXYZ(ex, ey, ez); }
   /**
    * @brief Set the magnetic field (T)
    */
   void SetBField(const XYZVector &bField) { fBField = bField; }
   void SetBField(double bx, double by, double bz) { fBField.SetXYZ(bx, by, bz); }

   /**
    * @brief Set the state of the particle.
    *
    * @param pos Position of the particle in mm.
    * @param mom Momentum of the particle in MeV/c.
    */
   void SetState(const XYZVector &pos, const XYZVector &mom)
   {
      fPos = pos;
      fMom = mom;
   }
   XYZVector GetPosition() const { return fPos; }
   XYZVector GetMomentum() const { return fMom; }

   /**
    * @brief Propagate the particle to the point of closest approach to the given point.
    *
    * Propagate to a given point in space, adjusting the magnitude of the stopping power
    * to ensure that a specific about of energy is lost during the propagation.
    *
    * @param point The point to approach.
    * @param eLoss If not 0, constrain the energy loss to this value (adjusting the stopping power).
    */
   void PropagateToPoint(const XYZVector &point, double eLoss = 0);

   /**
    * @brief Calculate the force acting on the particle.
    *
    * @param pos Position of the particle in mm.
    * @param mom Momentum of the particle in MeV/c.
    * @return The force acting on the particle in N.
    */
   XYZVector Force(XYZVector pos, XYZVector mom) const;

protected:
   /**
    * @brief Perform a single RK4 step for propagation.
    *
    * This method performs a single Runge-Kutta 4th order step to propagate the particle's state.
    * Updates fPos and fMom.
    */
   void RK4Step();
};

} // namespace AtTools
#endif // #ifndef ATPROPAGATOR_H
