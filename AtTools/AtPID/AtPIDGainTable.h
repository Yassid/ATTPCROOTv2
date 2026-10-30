#ifndef ATPIDGAINTABLE_H
#define ATPIDGAINTABLE_H

#include <istream>
#include <map>
#include <string>

namespace AtTools {

/// Per-run dE/dx gain factors for the PID estimators (AtSpyralPID / AtHelix3D::SetDedxGain).
///
/// WHY. The AT-TPC gain is not constant over a campaign: on a1975 D2 the deuteron band moves by
/// 0.12 in ln(sqrt(dE/dx) Brho) from run 0016 to 0097 (gain -20 %), seen identically by both
/// estimators. Pooled runs smear every PID band and a single gate is off-centre for some runs, a
/// run-dependent efficiency. Correcting per run narrowed the proton and deuteron bands by 12-13 %.
///
/// FORMAT (text): one line per run, `run gainSpyral gain3D`, factors multiply sqrt(dE/dx); `#` starts
/// a comment. Factors are relative to a reference run (1.0) and are produced from the data itself
/// (macro/Unpack_HDF5/a1975/D2_UKF/pid_gain_calib.C). A run that is NOT in the table returns -1:
/// a production must fail loudly rather than silently run uncorrected.
class AtPIDGainTable {
public:
   bool Load(const std::string &file);
   bool Load(std::istream &in);
   bool Has(int run) const { return fGain.count(run) > 0; }
   double GetSpyral(int run) const; ///< -1 if the run is absent
   double Get3D(int run) const;     ///< -1 if the run is absent
   std::size_t Size() const { return fGain.size(); }
   /// Run number from a file or run name such as "run_0016_multifit" (first digits after "run_");
   /// -1 if there is none.
   static int RunFromName(const std::string &name);

private:
   std::map<int, std::pair<double, double>> fGain;
};

} // namespace AtTools
#endif
