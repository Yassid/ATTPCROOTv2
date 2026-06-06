#ifndef ATPARTICLEID_H
#define ATPARTICLEID_H

#include "AtCut2D.h"

#include <Rtypes.h>

#include <string>

namespace AtTools {

/**
 * @brief A particle identification, C++ port of spyral_utils' ParticleID.
 *
 * Couples a 2D gate (AtCut2D, typically in the sqrt_dEdx vs brho plane) with a
 * nuclear species (Z, A). A track whose (x, y) observables fall inside the cut
 * is identified as this nucleus. The JSON read/written matches the spyral_utils
 * ParticleID format: the Cut2D fields plus top-level "Z" and "A":
 *
 *   {
 *       "name": "proton_cut",
 *       "xaxis": "sqrt_dEdx",
 *       "yaxis": "brho",
 *       "vertices": [[x1, y1], ...],
 *       "Z": 1,
 *       "A": 1
 *   }
 *
 * Not persisted to ROOT files — plain C++ types (per CLAUDE.md).
 */
class AtParticleID {
public:
   AtParticleID() = default;
   AtParticleID(AtCut2D cut, int Z, int A);

   /// True if (x, y) falls inside the gate.
   bool IsInside(double x, double y) const { return fCut.IsInside(x, y); }

   const AtCut2D &GetCut() const { return fCut; }
   int GetZ() const { return fZ; }
   int GetA() const { return fA; }
   double GetMassMeV() const { return fMassMeV; }   // nuclear mass, MeV/c^2
   double GetMassAmu() const;                        // nuclear mass, atomic mass units
   const std::string &GetName() const { return fCut.GetName(); }
   const std::string &GetXAxis() const { return fCut.GetXAxis(); }
   const std::string &GetYAxis() const { return fCut.GetYAxis(); }
   bool IsValid() const { return fCut.IsValid() && fA > 0; }

   static AtParticleID LoadJSON(const std::string &path);
   bool WriteJSON(const std::string &path) const;

   /// Nuclear mass (MeV/c^2) for a (Z, A). Exact for common light ions
   /// (p, d, t, 3He, 4He, 6/7Li, 12/13/14/15/16/17C); falls back to A*u for others.
   static double LookupMassMeV(int Z, int A);

private:
   AtCut2D fCut;
   int fZ{0};
   int fA{0};
   double fMassMeV{0.0};
};

} // namespace AtTools

#endif // ATPARTICLEID_H
