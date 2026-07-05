#ifndef ATBETHEBLOCHPID_H
#define ATBETHEBLOCHPID_H

#include <string>
#include <vector>

class AtTrack;

namespace AtTools {

/**
 * @brief Bethe-Bloch dE/dx vs rigidity particle identification.
 *
 * Classifies a pattern-recognized track among a set of mass hypotheses by comparing
 * its measured dE/dx to the Bethe-Bloch expectation at its magnetic RIGIDITY. This is
 * the standard TPC PID and is the right discriminant when a kinematic fit cannot tell
 * the species apart (e.g. K vs pi at a few hundred MeV/c over a short arc, where the
 * helices are near-identical): at fixed rigidity the heavier species sits higher on
 * the Bethe-Bloch curve and ionizes more.
 *
 * Momentum here is deliberately MASS-INDEPENDENT: the transverse rigidity
 * p_T = 0.299792458 * B[T] * R[mm] from the PRA circle radius, never a mass-dependent
 * fit KE (that would be circular). The dE/dx scale is detector/gain specific, so one
 * calibration constant k maps the Bethe-Bloch shape to the measured dE/dx units;
 * calibrate it from a known population (SetCalibration / SetReferenceCalibration).
 *
 * Not persisted to ROOT files — plain C++ types (per CLAUDE.md).
 */
class AtBetheBlochPID {
public:
   struct Species {
      std::string name; ///< base name, e.g. "pi", "K"
      double massMeV;   ///< rest mass, MeV/c^2
      int absPDG;       ///< |PDG| (sign applied from the track charge), e.g. 211, 321
   };

   AtBetheBlochPID() = default;

   /// Register a mass hypothesis. Order does not matter; Classify returns the best.
   void AddSpecies(const std::string &name, double massMeV, int absPDG);

   /// dE/dx-scale calibration: measured_dEdx = k * BetheBlochShape(betagamma).
   void SetCalibration(double k) { fK = k; }
   /// Convenience: derive k so that a reference species at reference rigidity/dE/dx matches.
   void SetReferenceCalibration(double refRigidityMeV, double refDeDx, double refMassMeV);
   double GetCalibration() const { return fK; }

   /// Mean-energy-loss Bethe-Bloch shape (arbitrary units; scale absorbed by k).
   static double BetheBlochShape(double p_MeV, double massMeV, double I_MeV = 40e-6);

   /// Expected measured dE/dx for a given rigidity and mass (= k * shape).
   double ExpectedDeDx(double rigidityMeV, double massMeV) const;

   /// Index of the best-matching registered species (nearest expected dE/dx in log
   /// space), or -1 if none registered / inputs invalid. betagamma uses rigidity/mass.
   int Classify(double rigidityMeV, double dEdx) const;

   const std::vector<Species> &GetSpecies() const { return fSpecies; }

   /// Track observables with the validated recipe: transverse rigidity from the PRA
   /// circle radius and the truncated-mean dE/dx (per-hit dQ/dl, drop the top
   /// @p trimHighFrac as the Landau tail). Returns false if the track is unusable.
   static bool Observables(const AtTrack &track, double bField_T, double &rigidityMeV, double &dEdx,
                           double trimHighFrac = 0.30);

private:
   std::vector<Species> fSpecies;
   double fK{1.0};
};

} // namespace AtTools

#endif // ATBETHEBLOCHPID_H
