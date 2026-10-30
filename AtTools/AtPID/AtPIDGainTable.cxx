#include "AtPIDGainTable.h"

#include <FairLogger.h>

#include <cctype>
#include <fstream>
#include <sstream>

using namespace AtTools;

bool AtPIDGainTable::Load(const std::string &file)
{
   std::ifstream in(file);
   if (!in) {
      LOG(error) << "AtPIDGainTable: cannot open " << file;
      return false;
   }
   return Load(in);
}

bool AtPIDGainTable::Load(std::istream &in)
{
   fGain.clear();
   std::string line;
   int nLine = 0;
   while (std::getline(in, line)) {
      ++nLine;
      auto hash = line.find('#');
      if (hash != std::string::npos)
         line.erase(hash);
      std::istringstream ss(line);
      int run;
      double gS, g3;
      if (!(ss >> run))
         continue; // blank or comment-only
      if (!(ss >> gS >> g3) || gS <= 0 || g3 <= 0) {
         LOG(error) << "AtPIDGainTable: malformed line " << nLine << ": '" << line << "'";
         fGain.clear();
         return false;
      }
      fGain[run] = {gS, g3};
   }
   return !fGain.empty();
}

double AtPIDGainTable::GetSpyral(int run) const
{
   auto it = fGain.find(run);
   return it == fGain.end() ? -1 : it->second.first;
}

double AtPIDGainTable::Get3D(int run) const
{
   auto it = fGain.find(run);
   return it == fGain.end() ? -1 : it->second.second;
}

int AtPIDGainTable::RunFromName(const std::string &name)
{
   auto p = name.rfind("run_");
   if (p == std::string::npos)
      return -1;
   p += 4;
   std::string digits;
   while (p < name.size() && std::isdigit(static_cast<unsigned char>(name[p])))
      digits += name[p++];
   return digits.empty() ? -1 : std::stoi(digits);
}
