#ifndef ATCUT2D_H
#define ATCUT2D_H

#include <Rtypes.h>

#include <string>
#include <utility>
#include <vector>

namespace AtTools {

/**
 * @brief A 2D polygon gate, C++ port of spyral_utils' Cut2D.
 *
 * Holds an ordered list of (x, y) vertices defining a closed polygon, plus the
 * names of the x/y observables it gates on (e.g. "sqrt_dEdx" and "brho"). The
 * JSON read/written matches the spyral_utils Cut2D format exactly, so gate files
 * produced by Spyral's particle_id notebook load directly:
 *
 *   {
 *       "name": "proton_cut",
 *       "xaxis": "sqrt_dEdx",
 *       "yaxis": "brho",
 *       "vertices": [[x1, y1], [x2, y2], ...]
 *   }
 *
 * Not persisted to ROOT files — plain C++ types (per CLAUDE.md).
 */
class AtCut2D {
public:
   AtCut2D() = default;
   AtCut2D(std::string name, std::vector<std::pair<double, double>> vertices, std::string xAxis = "DefaultAxis",
           std::string yAxis = "DefaultAxis");

   /// Point-in-polygon test (ray-casting). Mirrors spyral_utils Cut2D.is_inside.
   bool IsInside(double x, double y) const;

   const std::string &GetName() const { return fName; }
   const std::string &GetXAxis() const { return fXAxis; }
   const std::string &GetYAxis() const { return fYAxis; }
   const std::vector<std::pair<double, double>> &GetVertices() const { return fVertices; }
   bool IsValid() const { return fVertices.size() >= 3; }

   void SetName(std::string name) { fName = std::move(name); }
   void SetAxes(std::string xAxis, std::string yAxis)
   {
      fXAxis = std::move(xAxis);
      fYAxis = std::move(yAxis);
   }
   void SetVertices(std::vector<std::pair<double, double>> v) { fVertices = std::move(v); }

   /// Load a cut from a spyral_utils-format JSON file. Returns an invalid cut
   /// (IsValid()==false) on failure.
   static AtCut2D LoadJSON(const std::string &path);
   /// Write the cut to a spyral_utils-format JSON file. Returns false on failure.
   bool WriteJSON(const std::string &path) const;

private:
   std::string fName{"DefaultCut"};
   std::string fXAxis{"DefaultAxis"};
   std::string fYAxis{"DefaultAxis"};
   std::vector<std::pair<double, double>> fVertices;
};

} // namespace AtTools

#endif // ATCUT2D_H
