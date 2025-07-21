#ifndef ATPROPAGATOR_H
#define ATPROPAGATOR_H

#include "AtELossModel.h"

#include <functional>

#include "Math/Plane3D.h"
#include "Math/Point3D.h"
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
   using XYZPoint = ROOT::Math::XYZPoint;
   using Plane3D = ROOT::Math::Plane3D;
   using DistanceFunc = std::function<double(const XYZPoint &)>;
   XYZVector fEField{0, 0, 0}; // Electric field vector
   XYZVector fBField{0, 0, 0}; // Magnetic field vector

   const double fQ;                                 // Charge of the particle in Coulombs
   const double fMass;                              // Mass of the particle in MeV/c^2
   const std::unique_ptr<AtELossModel> fELossModel; // Energy loss model

   // Internal state variables for the propagator
   double fH = 1e-4;            /// Step size for propagation in m
   double fDelta = 1e-3;        /// Relative error tolerance for adaptive step size. 10^-3 means each 1m of propagation
                                /// introduces at most 1mm of error.
   double fETol = 1e-4;         /// Energy tolerance for convergence when fixing energy loss
   double fStopTol = 0.01;      /// Maximum kinetic energy to consider the particle stopped (MeV)
   double fDistTol = 1e-2;      /// Distance tolerance when considering positions equal. (mm)
   double fScalingFactor = 1.0; /// Scaling factor for energy loss

   XYZPoint fPos;  // Current position of the particle in mm
   XYZVector fMom; // Current momentum of the particle in MeV/c

   XYZPoint fLastPos;
   XYZVector fLastMom;

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
   void SetState(const XYZPoint &pos, const XYZVector &mom)
   {
      fPos = pos;
      fMom = mom;
   }

   void SetDelta(double delta) { fDelta = delta; }
   void SetH(double h) { fH = h; }

   XYZPoint GetPosition() const { return fPos; }
   XYZVector GetMomentum() const { return fMom; }

   /**
    * @brief Propagate the particle to the point of closest approach to the given point.
    *
    * Propagate to a given point in space, adjusting the magnitude of the stopping power
    * to ensure that a specific about of energy is lost during the propagation.
    *
    * @param point The point to approach.
    * @param eLoss If not 0, constrain the energy loss to this value by adjusting fScalingFactor.
    */
   void PropagateToPoint(const XYZPoint &point, double eLoss);

   /**
    * @brief Propagate the particle to the point of closest approach to the given point.
    *
    * @param point The point to approach.
    */
   void PropagateToPoint(const XYZPoint &point);
   void PropagateToPointAdaptive(const XYZPoint &point);

   /**
    * @brief Propagate the particle to the given plane.
    *
    * @param plane The plane to approach.
    */
   void PropagateToPlane(const Plane3D &plane);

   /**
    * @brief Calculate the force acting on the particle.
    *
    * @param pos Position of the particle in mm.
    * @param mom Momentum of the particle in MeV/c.
    * @return The force acting on the particle in N.
    */
   XYZVector Force(XYZPoint pos, XYZVector mom) const;

   /**
    * @brief Calculate the derivate of the momentum w.r.t. arc length.
    *
    * @param pos Position of the particle in mm.
    * @param mom Momentum of the particle in MeV/c.
    * @return The derivative of the momentum w.r.t. arc length in N/m.
    */
   XYZVector dpds(const XYZPoint &pos, const XYZVector &mom) const;

   /**
    * @brief Calculate the second derivative of the position w.r.t. arc length.
    *
    * \frac{d^2\vec{x}}{ds^2} = \frac{1}{p} \left( \frac{d\vec{p}}{ds} - \hat{p} (\hat{p} \cdot \frac{d\vec{p}}{ds})
    * \right)
    *
    * @param pos Position of the particle in mm.
    * @param mom Momentum of the particle in MeV/c.
    * @return The second derivative of the position w.r.t. arc length in m/m^2.
    */
   XYZVector d2xds2(const XYZPoint &pos, const XYZVector &mom) const;

   XYZVector dxds(const XYZPoint &pos, const XYZVector &mom) const
   {
      return mom.Unit(); // The derivative of the position is just the unit vector of the momentum.
   }

protected:
   /**
    * @brief Perform a single RK4 step for propagation.
    *
    * This method performs a single Runge-Kutta 4th order step to propagate the particle's state.
    * Updates fPos and fMom.
    * @param h Step size for the RK4 step in meters.
    */
   void RK4Step(double h);

   /**
    * @brief Perform a single RK4 step using the Nystrom method.
    * This method performs a single Runge-Kutta 4th order step using the Nystrom method
    * to propagate the particle's state. Updates fPos and fMom.
    * @param h Step size for the RK4 step in meters.
    */
   void RK4StepNystrom(double h);

   /**
    * @brief Perform an adaptive RK4 step for propagation.
    * This method performs an adaptive Runge-Kutta 4th order step to propagate the particle's state.
    * Updates fPos and fMom.
    *
    * The error is based on the difference in the positions at the end of the step. i.e:
    * \eps = \sqrt{\eps_x^2 + \eps_y^2 + \eps_z^2}, where
    * \eps_x = 1/30*|x_1 - x_2|, where x_1 is using h/2 and x_2 is using h.
    *
    * Step size is adjust to ensure the local error is less than fDelta.
    *
    * @param h Step size for the RK4 step in seconds. Modified in place to reflect the new step size.
    * @return True if the step was accepted, false otherwise.
    */
   bool RK4StepAdaptive(double &h);

   void PropagateTo(DistanceFunc distanceFunc);

   bool ReachedPOCA(const XYZPoint &point);
   bool IntersectedPlane(const Plane3D &plane);
};

static constexpr double a21 = 1.0 / 5.0;
static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                        a65 = -5103.0 / 18656.0;
static constexpr double a71 = 35.0 / 384.0, a72 = 0.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                        a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;
// b (5th-order)
static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0,
                        b6 = 11.0 / 84.0;
// b* (4th-order, “star”)
static constexpr double bs1 = 5179.0 / 57600.0, bs3 = 7571.0 / 16695.0, bs4 = 393.0 / 640.0, bs5 = -92097.0 / 339200.0,
                        bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

} // namespace AtTools
#endif // #ifndef ATPROPAGATOR_H
