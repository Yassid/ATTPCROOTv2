/// @file gain_C15p.h
/// @brief Load a C15p gain-match table and apply it at read time.
///
/// The reco persists RAW dE/dx on purpose, so the gain correction is applied where the plane is
/// consumed rather than baked into 72 GB of reco. Anything that reads the PID caches and shows or
/// cuts on dE/dx must apply the SAME table, or a gate drawn on a matched plane gets applied to
/// unmatched values and quietly selects the wrong tracks.
///
///     dEdx *= f(run)      sqrtdEdx *= sqrt(f(run))
///
/// Same table and same convention as AtGainMatchTask, which is what applies it inside the
/// framework chain when a production is run with gain matching on.

#ifndef GAIN_C15D_H
#define GAIN_C15D_H

#include <fstream>
#include <map>
#include <sstream>
#include <string>

/// run -> factor. Empty map = no matching, and callers must SAY SO rather than silently
/// presenting a raw plane as a matched one.
inline std::map<int, double> LoadGainTable_C15p(const TString &path, bool verbose = true)
{
   std::map<int, double> table;
   if (path.Length() == 0)
      return table;
   std::ifstream in(path.Data());
   if (!in) {
      std::cout << "\033[1;31mERROR: gain table " << path << " not found -- plane stays RAW.\033[0m\n";
      return table;
   }
   std::string line;
   int nMeas = 0, nInterp = 0;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line.rfind("run,", 0) == 0)
         continue;
      std::stringstream ss(line);
      std::string r, f, src;
      if (!std::getline(ss, r, ',') || !std::getline(ss, f, ','))
         continue;
      std::getline(ss, src, ',');
      const double fac = std::atof(f.c_str());
      if (!(fac > 0))
         continue;
      table[std::atoi(r.c_str())] = fac;
      if (src == "measured")
         ++nMeas;
      else
         ++nInterp;
   }
   if (verbose)
      std::cout << "  gain table  : " << path << "  (" << table.size() << " runs: " << nMeas
                << " measured, " << nInterp << " interpolated/held)\n";
   return table;
}

/// Factor for one run. A run ABSENT from the table returns 1.0 and sets `missing` -- callers are
/// expected to count those and report them, because an unmatched run sitting inside a matched
/// sample is exactly the thing that is invisible downstream.
inline double GainFactor_C15p(const std::map<int, double> &table, int run, bool &missing)
{
   if (table.empty()) {
      missing = false; // no table at all is a deliberate raw plane, not a missing run
      return 1.0;
   }
   auto it = table.find(run);
   if (it == table.end()) {
      missing = true;
      return 1.0;
   }
   missing = false;
   return it->second;
}

#endif // GAIN_C15D_H
