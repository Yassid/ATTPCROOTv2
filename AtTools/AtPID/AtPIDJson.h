#ifndef ATPIDJSON_H
#define ATPIDJSON_H

// Minimal JSON extraction helpers shared by AtCut2D and AtParticleID.
// Header-only, inline, NOT registered in the ROOT dictionary. Scoped to the
// fixed spyral_utils Cut2D / ParticleID schema (flat object with string keys,
// integer keys, and one "vertices" array of [x, y] pairs) — not a general JSON
// parser, but robust to whitespace, ordering and scientific notation in the
// machine-generated gate files.

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace AtTools {
namespace PIDJson {

/// Read an entire file into a string. Empty string on failure.
inline std::string ReadFile(const std::string &path)
{
   std::ifstream in(path);
   if (!in)
      return "";
   std::stringstream ss;
   ss << in.rdbuf();
   return ss.str();
}

/// Value of a "key": "value" string field. Empty if not found.
inline std::string ExtractString(const std::string &content, const std::string &key, const std::string &def = "")
{
   const std::string needle = "\"" + key + "\"";
   auto k = content.find(needle);
   if (k == std::string::npos)
      return def;
   auto colon = content.find(':', k + needle.size());
   if (colon == std::string::npos)
      return def;
   auto q1 = content.find('"', colon + 1);
   if (q1 == std::string::npos)
      return def;
   auto q2 = content.find('"', q1 + 1);
   if (q2 == std::string::npos)
      return def;
   return content.substr(q1 + 1, q2 - q1 - 1);
}

/// Value of a "key": <number> integer field. `def` if not found.
inline int ExtractInt(const std::string &content, const std::string &key, int def = 0)
{
   const std::string needle = "\"" + key + "\"";
   auto k = content.find(needle);
   if (k == std::string::npos)
      return def;
   auto colon = content.find(':', k + needle.size());
   if (colon == std::string::npos)
      return def;
   return static_cast<int>(std::strtol(content.c_str() + colon + 1, nullptr, 10));
}

/// All [x, y] pairs inside the "vertices" array.
inline std::vector<std::pair<double, double>> ExtractVertices(const std::string &content)
{
   std::vector<std::pair<double, double>> verts;
   auto k = content.find("\"vertices\"");
   if (k == std::string::npos)
      return verts;
   auto open = content.find('[', k);
   if (open == std::string::npos)
      return verts;
   // Find the matching close bracket of the outer array (depth counting).
   int depth = 0;
   std::size_t close = open;
   for (std::size_t i = open; i < content.size(); ++i) {
      if (content[i] == '[')
         ++depth;
      else if (content[i] == ']') {
         --depth;
         if (depth == 0) {
            close = i;
            break;
         }
      }
   }
   std::string body = content.substr(open, close - open + 1);
   // Replace structural characters with spaces, then stream all doubles in order.
   for (char &c : body)
      if (c == '[' || c == ']' || c == ',')
         c = ' ';
   std::stringstream ss(body);
   std::vector<double> nums;
   double v;
   while (ss >> v)
      nums.push_back(v);
   for (std::size_t i = 0; i + 1 < nums.size(); i += 2)
      verts.emplace_back(nums[i], nums[i + 1]);
   return verts;
}

} // namespace PIDJson
} // namespace AtTools

#endif // ATPIDJSON_H
