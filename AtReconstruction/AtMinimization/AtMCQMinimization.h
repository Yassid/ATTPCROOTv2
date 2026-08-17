/*******************************************************************
 * Monte Carlo minimization using the charge deposited on the pads  *
 * Original implementation (FairRootv18.00 branch): 29-07-2016      *
 * Author: Y. Ayyad and W. Mittig (NSCL)                            *
 ********************************************************************/

#ifndef ATMCQMINIMIZATION_H
#define ATMCQMINIMIZATION_H

#include "AtMinimization.h"

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>

#include <functional>
#include <memory>
#include <vector>

class AtDigiPar;
class AtHit;
class AtMap;

namespace AtTools {
class AtELossModel;
}

namespace MCMinimization {

/**
 * @brief Monte Carlo minimizer comparing the charge collected by each pad.
 *
 * The track described by the seed is propagated through the gas (Euler integration of the
 * Lorentz force plus the energy loss), the primary electrons are spread over the pad plane
 * with a simple transverse straggling model, and the resulting charge pattern is compared
 * with the experimental one. The comparison drives a random walk with a shrinking step
 * (see VaryState()) that keeps every trial improving the objective function.
 *
 * The objective function is the average of
 *  - a charge term, comparing the charge and the charge weighted z of every pad (Chi2Q()), and
 *  - either a position term, comparing the center of gravity of the hits of each time bucket
 *    with the simulated trajectory (Chi2Pos(), enabled by default), or a range term comparing
 *    the end point of the track with the simulated one (Chi2Range(), SetUseRangeChi2()).
 *
 * Minimal setup, taking the energy loss from CATIMA:
 * @code
 *    auto eLoss = std::make_shared<AtTools::AtELossCATIMA>(density, material);
 *    eLoss->SetProjectile(A, Z, massAmu);
 *
 *    auto min = MCMinimization::AtMCQMinimization();
 *    min.SetMap(fMap);          // pad plane geometry, mandatory
 *    min.SetParticle(A, Z);     // fitted particle
 *    min.SetELossModel(eLoss);  // stopping power and range to energy conversion
 *    min.Minimize(seed, event);
 * @endcode
 *
 * The charge response (ionization energy, transverse diffusion, micromegas and GET gain) is read
 * from AtDigiPar, so the track simulated here deposits the same charge on the same pads as a
 * track pushed through AtClusterize and AtPulse. Only the overall scale can be overriden:
 * @code
 *    min.SetGain(adcPerElectron);
 *    auto res = min.GetFitPar();
 * @endcode
 *
 * Any other AtTools::AtELossModel (a table, SRIM, ...) works the same way. The parametrized
 * dE/dx and range to energy functions of the FairRootv18.00 macros are still accepted, see
 * SetStoppingPower() and SetRangeToEnergy().
 *
 * Differences with the FairRootv18.00 implementation are listed at the top of the source file.
 */
class AtMCQMinimization : public AtMinimization {
public:
   using MapPtr = std::shared_ptr<AtMap>;
   using ELossModelPtr = std::shared_ptr<AtTools::AtELossModel>;

   /// Stopping power dE/dx [MeV/cm] in the gas of the parameter file, as a function of the
   /// kinetic energy [MeV]. The random walk scales it by the (dimensionless) density of MCState.
   using StoppingPowerFunc = std::function<double(double)>;
   /// Kinetic energy [MeV] as a function of the range [mm] in the gas of the detector.
   using RangeToEnergyFunc = std::function<double(double)>;
   /// Signature of the energy loss functions of the macros of the FairRootv18.00 branch.
   using ParametrisedFunc = std::function<double(double, std::vector<double> &)>;

   /**
    * @brief The parameters explored by the random walk.
    *
    * The position is in the pad plane frame and in cm (the seed is in mm), the angles are in
    * radians and the magnetic field is in gauss, as in the original implementation.
    */
   struct MCState {
      double fX{0};            //< [cm]
      double fY{0};            //< [cm]
      double fZ{0};            //< [cm]
      double fTheta{0};        //< [rad]
      double fPhi{0};          //< [rad]
      double fB{0};            //< [G]
      double fDensityScale{1}; //< Scaling of the stopping power, 1 = the gas of the parameter file
      double fRadius{0};       //< Radius of curvature [mm], or the range [mm] when the field is off
   };

private:
   /// Number of Euler steps of the simulation, and size of the per-step buffers.
   static constexpr int kMaxIntegrationSteps = 10000;
   /// Sentinel marking a time bucket the simulated track never reached.
   static constexpr double kNoPoint = -10000;

   MapPtr fMap{}; //< Pad plane geometry, mandatory

   // Parameters read from AtDigiPar by Init()
   double fDriftVelocity{0}; //< [cm/us]
   int fTBTime{0};           //< Time bucket length [ns]
   double fBField{0};        //< [T]
   double fZPadPlane{0};     //< Position of the pad plane along the beam axis [mm]
   double fDensity{0};       //< Gas density, in the units of the parameter file
   double fPressure{0};      //< Gas pressure [torr]
   int fEntTB{0};            //< Time bucket of the entrance window

   /* Charge response, the same parameters the digitization uses (AtClusterize and AtPulse), so
    * that a track simulated here and a track pushed through the digitization chain deposit the
    * same charge on the same pads. */
   double fEIonize{0};   //< Effective ionization energy of the gas [MeV]
   double fCoefT{0};     //< Transverse diffusion coefficient [cm2/us]
   double fCoefL{0};     //< Longitudinal diffusion coefficient [cm2/us], see SimulateTrack()
   double fGain{0};      //< ADC per primary electron, micromegas gain times GET gain
   double fMaxRange{0};  //< Length of the active volume along the beam axis [mm]
   bool fGainSet{false}; //< True once SetGain() was called, so Init() does not overwrite it

protected:
   /* Geometry of the tilted detector of the 2016 campaigns: the pad plane rotation, the tilt of
    * the detector inside the solenoid and the Lorentz angle of the drift. The current AtDigiPar
    * does not carry them anymore and the rest of the reconstruction assumes they are zero, so
    * they are not exposed; a derived class working with that geometry can still set them and
    * every frame transformation of the class honours them. */
   double fTiltAng{0};      //< Tilt of the detector with respect to the beam axis [rad]
   double fThetaLorentz{0}; //< Lorentz angle of the drifting electrons [rad]
   double fThetaRot{0};     //< Azimuthal orientation of the Lorentz drift [rad]
   double fThetaPad{0};     //< Rotation of the pad plane [rad]

private:
   // Particle
   int fA{0};               //< Mass number
   int fZ{0};               //< Atomic number
   double fRestMass{0};     //< [MeV]
   double fChargeToMass{0}; //< Charge over mass, in the units of the equation of motion

   ELossModelPtr fELossModel{};        //< Set by SetELossModel(), owns what the two functions below use
   StoppingPowerFunc fStoppingPower{}; //< dE/dx [MeV/cm] vs kinetic energy [MeV]
   RangeToEnergyFunc fRangeToEnergy{}; //< Kinetic energy [MeV] vs range [mm] in the gas

   // Search configuration
   std::vector<double> fStepPar{4.0, 4.0, 0.1, 0.5, 0.5, 0.5, 0.0, 0.0}; //< Default steps (46Ar)
   /* Widths entering the objective function. They are not resolutions: they set how much a
    * disagreement in charge weighs against a disagreement in position, so they are the handle
    * on the balance between the two terms. */
   double fSigmaQ{0.2}; //< Charge width, as a fraction of the sum Qsim + Qexp
   double fSigmaZ{4.0}; //< Width of the charge weighted position of a pad [mm]

   int fNumCoarseIter{5};     //< Number of step reductions
   int fNumFineIter{40};      //< Number of trials per step
   int fIntegrationSteps{10}; //< Euler steps per time bucket

   bool fUsePosChi2{true};    //< Add the position term to the objective function
   bool fUseRangeChi2{false}; //< Replace the position term by the range term
   bool fUseGeoVertex{false}; //< Seed z from the vertex position instead of from its time bucket
   bool fBackwardProp{true};  //< Extrapolate the track back to the beam axis to find the vertex
   bool fVerbose{true};

   double fEntZ0{0};      //< Position of the entrance window (cathode) [mm], used with fUseGeoVertex
   bool fEntTBSet{false}; //< True once SetEntTB() was called, so Init() does not overwrite it

   // State of the current fit
   bool fIsInit{false};
   double fDzStep{0};      //< Length of a time bucket along the drift direction [cm]
   double fBFieldGauss{0}; //< fBField in gauss, the working unit of the equation of motion
   double fBeamRange{0};   //< Range of the beam particle, used when the vertex is not extrapolated
   int fChi2Points{0};     //< Number of pads entering the charge chi2
   int fLastTB{0};         //< Last time bucket reached by the last simulation

   /* Experimental track, filled by Minimize() */
   std::vector<double> fQExp;                         //< Experimental charge per pad
   std::vector<double> fZExp;                         //< z of the experimental hit of each pad [mm]
   std::vector<XYZPoint> fExpTrack;                   //< Position of every experimental hit [mm]
   std::vector<int> fExpTB;                           //< Time bucket of every experimental hit
   std::vector<std::vector<const AtHit *>> fHitsByTB; //< Hits sorted by descending time bucket

   /* Simulated track, filled by SimulateTrack() */
   std::vector<double> fQSim;        //< Simulated charge per pad
   std::vector<double> fZSim;        //< Charge weighted z of each pad [mm]
   std::vector<XYZPoint> fSimTrack;  //< Simulated trajectory, one point per integration step [mm]
   std::vector<XYZPoint> fSimPads;   //< Center of the pads above threshold, z is the centroid [mm]
   std::vector<double> fSimCharge;   //< Charge of the pads of fSimPads, before the gain
   std::vector<XYZPoint> fTBTrack;   //< Simulated trajectory indexed by time bucket [mm]
   std::vector<XYZPoint> fBackTrack; //< Backward extrapolation of the fitted track [mm]
   double fSimEnergy{0};             //< Energy of the last simulated track [MeV/u]

   /* Scratch buffers of SimulateTrack(), members to keep them out of the hot loop */
   std::vector<double> fStepX;
   std::vector<double> fStepY;
   std::vector<double> fStepZ;
   std::vector<double> fStepQ;

public:
   AtMCQMinimization() = default;
   ~AtMCQMinimization() override = default;

   /// Read AtDigiPar and size the internal buffers. Called by Minimize() if not called before.
   bool Init();

   void SetMap(MapPtr map) { fMap = std::move(map); }
   void SetParticle(int A, int Z);

   /**
    * @brief Take the energy loss from one of the models of AtTools, e.g. AtELossCATIMA.
    *
    * Sets both the stopping power and, by inverting the range of the model, the range to energy
    * conversion the fit needs when the magnetic field is off. The projectile of the model has to
    * be the particle passed to SetParticle(), and the model has to be built for the gas of the
    * detector: its density is what fixes the absolute energy loss, the density of the parameter
    * file is only used as the starting point of the (by default frozen) density scan.
    */
   void SetELossModel(ELossModelPtr model);

   void SetStoppingPower(StoppingPowerFunc func) { fStoppingPower = std::move(func); }
   void SetRangeToEnergy(RangeToEnergyFunc func) { fRangeToEnergy = std::move(func); }
   /**
    * @brief Overload taking the (energy, parameters) dE/dx functions of the FairRootv18.00 macros.
    *
    * Those return the mass stopping power, so the result is multiplied by the density of the
    * parameter file, as the original implementation did.
    */
   void SetStoppingPower(ParametrisedFunc func, std::vector<double> par);
   /**
    * @brief Overload taking the (range, parameters) functions of the FairRootv18.00 macros.
    *
    * Those are defined at 760 torr, so the range is scaled by the pressure of the gas before
    * being passed to the function, as the original implementation did.
    */
   void SetRangeToEnergy(ParametrisedFunc func, std::vector<double> par);

   /**
    * @brief Steps of the random walk, i.e. how far a trial track can move from the current one.
    *
    * In order: theta [deg], phi [deg], radius [relative], x, y, z [cm], magnetic field and
    * density [relative]. Every step is scaled down by 1.4 at each of the coarse iterations.
    * The defaults are the ones the original implementation used for 46Ar; a vector longer than
    * eight entries is accepted, so the ten entry arrays of the old macros still work, and the
    * extra entries are ignored.
    */
   void SetStepParameters(const std::vector<double> &par);
   void SetNumIterations(int coarse, int fine);
   void SetIntegrationSteps(int steps) { fIntegrationSteps = steps; }

   /**
    * @brief Override the ADC per primary electron of the parameter file.
    *
    * By default the charge of the simulated track is scaled the way AtPulse scales it, i.e. by
    * the micromegas gain times the gain of the GET electronics. Since the objective function
    * compares charges through (Qsim - Qexp) / (Qsim + Qexp) it is sensitive to that scale, so
    * a run whose amplitudes are not reproduced by the nominal gain can set it here.
    */
   void SetGain(double adcPerElectron)
   {
      fGain = adcPerElectron;
      fGainSet = true;
   }
   double GetGain() const { return fGain; }

   void SetEntTB(int value)
   {
      fEntTB = value;
      fEntTBSet = true;
   }
   void SetEntZ0(double value) { fEntZ0 = value; }
   void SetUseGeoVertex(bool value = true) { fUseGeoVertex = value; }
   void SetBackwardPropagation(bool value = true) { fBackwardProp = value; }
   void SetUsePosChi2(bool value = true) { fUsePosChi2 = value; }
   void SetUseRangeChi2(bool value = true) { fUseRangeChi2 = value; }
   void SetVerbose(bool value = true) { fVerbose = value; }

   /// Length of the active volume [mm], only used when the vertex is not extrapolated back.
   void SetMaxRange(double mm) { fMaxRange = mm; }

   /// Charge width of the objective function, as a fraction of the sum Qsim + Qexp.
   void SetSigmaQ(double value) { fSigmaQ = value; }
   /// Width of the charge weighted position of a pad [mm] in the objective function.
   void SetSigmaZ(double value) { fSigmaZ = value; }
   double GetSigmaQ() const { return fSigmaQ; }
   double GetSigmaZ() const { return fSigmaZ; }

   using AtMinimization::Minimize;
   bool Minimize(const TrackSeed &seed, const HitVector &hits) override;
   void ResetParameters() override;

   /// Trajectory of the fitted track [mm], one point per integration step (SetIntegrationSteps()).
   const std::vector<XYZPoint> &GetSimTrack() const { return fSimTrack; }
   /// Center of the pads fired by the fitted track, z holds the charge weighted position [mm].
   const std::vector<XYZPoint> &GetSimPads() const { return fSimPads; }
   /// Charge of the pads returned by GetSimPads(), before the gain is applied.
   const std::vector<double> &GetSimCharge() const { return fSimCharge; }
   /// Position of the experimental hits that were fitted [mm].
   const std::vector<XYZPoint> &GetExpTrack() const { return fExpTrack; }
   /// Time bucket of the experimental hits that were fitted.
   const std::vector<int> &GetExpTB() const { return fExpTB; }
   /// Backward extrapolation of the fitted track towards the beam axis [mm].
   const std::vector<XYZPoint> &GetBackTrack() const { return fBackTrack; }

   /// Pad plane frame [cm] to laboratory frame [cm].
   XYZPoint TransformIniPos(const XYZPoint &pos) const;
   /// Laboratory frame [cm] to pad plane frame [cm].
   XYZPoint InvTransIniPos(const XYZPoint &pos) const;

protected:
   /// Seed the random walk, i.e. translate the seed into the frame of the simulation.
   MCState SeedState(const TrackSeed &seed) const;
   /// Draw a new state around the current one. The steps shrink with the iteration number.
   MCState VaryState(const MCState &state, int iteration) const;

   /**
    * @brief Simulate the track described by the state and collect its charge on the pad plane.
    *
    * Fills fQSim, fZSim, fSimTrack, fSimPads, fSimCharge, fTBTrack and fSimEnergy.
    */
   void SimulateTrack(const MCState &state);
   /// Extrapolate the fitted track back to the beam axis and fill the vertex of the result.
   void BackwardExtrapolation(const MCState &state);

   /// Charge and charge weighted position term of the objective function.
   /// Also updates the number of pads used, needed by Chi2Range().
   double Chi2Q();
   /// Position term, comparing the center of gravity of each time bucket with the simulation.
   double Chi2Pos(int iteration, int numExpPoints) const;
   /// Range term, comparing the end point of the experimental and simulated tracks.
   double Chi2Range(const HitVector &hits) const;

   /// ADC per primary electron, the conversion the digitization applies (AtPulse).
   static double GainFromParameters(const AtDigiPar &par);
   /// Kinetic energy per nucleon [MeV/u] of a particle of magnetic rigidity bRho [Tm].
   static double GetEnergy(double A, double Z, double bRho);
   /// Kinetic energy [MeV] of a particle of the given range [mm], by inverting fELossModel.
   double RangeToEnergyFromModel(double range) const;

   void ClearTrack();
   void FillExperimentalTrack(const HitVector &hits, int vertexTB);
   void PrintResult() const;
};

} // namespace MCMinimization

#endif // ATMCQMINIMIZATION_H
