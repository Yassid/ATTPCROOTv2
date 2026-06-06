#include "AtCut2D.h"

#include "AtPIDJson.h"

#include <FairLogger.h>

#include <fstream>
#include <iomanip>
#include <utility>

using namespace AtTools;

AtCut2D::AtCut2D(std::string name, std::vector<std::pair<double, double>> vertices, std::string xAxis,
                 std::string yAxis)
   : fName(std::move(name)), fXAxis(std::move(xAxis)), fYAxis(std::move(yAxis)), fVertices(std::move(vertices))
{
}

bool AtCut2D::IsInside(double x, double y) const
{
   // Standard ray-casting point-in-polygon. Mirrors the result of
   // shapely Polygon.contains used by spyral_utils Cut2D.
   const std::size_t n = fVertices.size();
   if (n < 3)
      return false;
   bool inside = false;
   for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
      const double xi = fVertices[i].first, yi = fVertices[i].second;
      const double xj = fVertices[j].first, yj = fVertices[j].second;
      const bool intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
      if (intersect)
         inside = !inside;
   }
   return inside;
}

AtCut2D AtCut2D::LoadJSON(const std::string &path)
{
   const std::string content = PIDJson::ReadFile(path);
   if (content.empty()) {
      LOG(error) << "AtCut2D::LoadJSON: could not read " << path;
      return {};
   }
   std::string name = PIDJson::ExtractString(content, "name", "DefaultCut");
   std::string xaxis = PIDJson::ExtractString(content, "xaxis", "DefaultAxis");
   std::string yaxis = PIDJson::ExtractString(content, "yaxis", "DefaultAxis");
   auto verts = PIDJson::ExtractVertices(content);
   if (verts.size() < 3) {
      LOG(error) << "AtCut2D::LoadJSON: fewer than 3 vertices parsed from " << path;
      return {};
   }
   LOG(info) << "AtCut2D::LoadJSON: loaded '" << name << "' (" << verts.size() << " vertices, " << xaxis << " vs "
             << yaxis << ") from " << path;
   return AtCut2D(std::move(name), std::move(verts), std::move(xaxis), std::move(yaxis));
}

bool AtCut2D::WriteJSON(const std::string &path) const
{
   std::ofstream out(path);
   if (!out) {
      LOG(error) << "AtCut2D::WriteJSON: cannot open " << path;
      return false;
   }
   out << std::setprecision(15);
   out << "{\n";
   out << "    \"name\": \"" << fName << "\",\n";
   out << "    \"xaxis\": \"" << fXAxis << "\",\n";
   out << "    \"yaxis\": \"" << fYAxis << "\",\n";
   out << "    \"vertices\": [\n";
   for (std::size_t i = 0; i < fVertices.size(); ++i) {
      out << "        [\n            " << fVertices[i].first << ",\n            " << fVertices[i].second << "\n        ]";
      out << (i + 1 < fVertices.size() ? ",\n" : "\n");
   }
   out << "    ]\n";
   out << "}\n";
   return true;
}
