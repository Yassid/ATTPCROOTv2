#include "AtParticleID.h"

#include "AtCut2D.h"
#include "AtPIDJson.h"

#include <FairLogger.h>

#include <fstream>
#include <iomanip>
#include <map>
#include <utility>

using namespace AtTools;

namespace {
constexpr double kU_MeV = 931.49410242; // 1 atomic mass unit in MeV/c^2
}

AtParticleID::AtParticleID(AtCut2D cut, int Z, int A)
   : fCut(std::move(cut)), fZ(Z), fA(A), fMassMeV(LookupMassMeV(Z, A))
{
}

double AtParticleID::GetMassAmu() const
{
   return fMassMeV / kU_MeV;
}

double AtParticleID::LookupMassMeV(int Z, int A)
{
   // Nuclear (not atomic) masses in MeV/c^2 for the common AT-TPC ejectiles and
   // carbon-beam isotopes. Keyed by (Z, A). Fallback: A * u (good to ~1%).
   static const std::map<std::pair<int, int>, double> table = {
      {{1, 1}, 938.27208816},  {{1, 2}, 1875.61294257},  {{1, 3}, 2808.92113298}, // p, d, t
      {{2, 3}, 2808.39160743}, {{2, 4}, 3727.3794066},                            // 3He, 4He
      {{3, 6}, 5601.518},      {{3, 7}, 6533.833},                                // 6Li, 7Li
      {{6, 12}, 11174.8625},   {{6, 13}, 12109.481},     {{6, 14}, 13040.870},    // 12C, 13C, 14C
      {{6, 15}, 13983.586},    {{6, 16}, 14914.529},     {{6, 17}, 15846.595},    // 15C, 16C, 17C
   };
   auto it = table.find({Z, A});
   if (it != table.end())
      return it->second;
   return A * kU_MeV;
}

AtParticleID AtParticleID::LoadJSON(const std::string &path)
{
   AtCut2D cut = AtCut2D::LoadJSON(path);
   if (!cut.IsValid()) {
      LOG(error) << "AtParticleID::LoadJSON: invalid cut in " << path;
      return {};
   }
   const std::string content = PIDJson::ReadFile(path);
   int Z = PIDJson::ExtractInt(content, "Z", 0);
   int A = PIDJson::ExtractInt(content, "A", 0);
   if (A <= 0) {
      LOG(warning) << "AtParticleID::LoadJSON: no valid (Z,A) in " << path
                   << " — loaded cut only, nucleus undefined.";
   }
   LOG(info) << "AtParticleID::LoadJSON: '" << cut.GetName() << "' -> Z=" << Z << " A=" << A
             << " (mass " << LookupMassMeV(Z, A) << " MeV)";
   return AtParticleID(std::move(cut), Z, A);
}

bool AtParticleID::WriteJSON(const std::string &path) const
{
   // Write the Cut2D body, then inject "Z" and "A" before the closing brace so
   // the file round-trips through spyral_utils' deserialize_particle_id.
   std::ofstream out(path);
   if (!out) {
      LOG(error) << "AtParticleID::WriteJSON: cannot open " << path;
      return false;
   }
   const auto &verts = fCut.GetVertices();
   out << std::setprecision(15);
   out << "{\n";
   out << "    \"name\": \"" << fCut.GetName() << "\",\n";
   out << "    \"xaxis\": \"" << fCut.GetXAxis() << "\",\n";
   out << "    \"yaxis\": \"" << fCut.GetYAxis() << "\",\n";
   out << "    \"vertices\": [\n";
   for (std::size_t i = 0; i < verts.size(); ++i) {
      out << "        [\n            " << verts[i].first << ",\n            " << verts[i].second << "\n        ]";
      out << (i + 1 < verts.size() ? ",\n" : "\n");
   }
   out << "    ],\n";
   out << "    \"Z\": " << fZ << ",\n";
   out << "    \"A\": " << fA << "\n";
   out << "}\n";
   return true;
}
