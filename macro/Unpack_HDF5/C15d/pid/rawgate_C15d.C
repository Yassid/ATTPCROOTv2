/// @file rawgate_C15d.C
/// @brief Re-express a gain-matched PID gate in RAW dE/dx for one run, so the FITTER can apply it.
///
///   root -b -q 'pid/rawgate_C15d.C("pid/deuteron_C15d.json", 26)'   -> pid/raw/deuteron_C15d_r26.json
///
/// WHY THIS EXISTS. Gating before the fit rather than after is worth roughly an order of magnitude
/// of CPU: on the C15d (d,d') pass, fitting every track over 75 runs meant ~1.9M fits to keep the
/// ~133k inside the deuteron gate. But `AtGenfitter::SetPIDGate` evaluates its own AtSpyralPID on
/// RAW dE/dx, while every gate in this workspace is drawn on the GAIN-MATCHED plane -- so handing
/// the matched gate straight to the fitter selects a different set of tracks, silently.
///
/// The gain match is a pure per-run scale, dEdx -> f(run) * dEdx, so on the sqrt axis the gate
/// transforms exactly:
///
///     x_raw = x_matched / sqrt(f(run))
///
/// Brho is untouched. So a per-run raw gate selects precisely the tracks the matched gate does,
/// with no approximation -- provided the fitter's AtSpyralPID computes dE/dx the same way
/// AtPIDTask did. That is an assumption about two instances of the same class, not a theorem:
/// VALIDATE IT against the (run, event, trackID) join before trusting a production to it
/// (pid/validate_rawgate_C15d.C).

#include "../gain_C15d.h"

void rawgate_C15d(TString gateFile = "pid/deuteron_C15d.json", Int_t run = 26,
                  TString gainTable = "gainmatch_C15d.csv", TString outDir = "pid/raw/")
{
   gSystem->Load("libAtTools.so");

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString gt = gainTable;
   if (gt.Length() && !gt.BeginsWith("/") && gSystem->AccessPathName(gt)) {
      gt = here + "/../" + gainTable;
      if (gSystem->AccessPathName(gt))
         gt = here + "/" + gainTable;
   }
   auto gain = LoadGainTable_C15d(gt, false);
   if (gain.empty()) {
      std::cout << "\033[1;31mERROR: no gain table -- cannot build a raw gate.\033[0m\n";
      return;
   }
   bool missing = false;
   const double f = GainFactor_C15d(gain, run, missing);
   if (missing) {
      std::cout << "\033[1;31mERROR: run " << run << " has no gain entry. A raw gate for it would be "
                << "the matched gate unscaled, i.e. the wrong tracks.\033[0m\n";
      return;
   }

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   const auto &cut = pid.GetCut();
   if (!cut.IsValid()) {
      std::cout << "\033[1;31mERROR: cannot load " << gateFile << "\033[0m\n";
      return;
   }

   const double s = 1.0 / std::sqrt(f);
   std::vector<std::pair<double, double>> v;
   for (const auto &p : cut.GetVertices())
      v.emplace_back(p.first * s, p.second); // Brho is unaffected by a dE/dx scale

   TString base = gSystem->BaseName(gateFile);
   base.ReplaceAll(".json", "");
   gSystem->mkdir(outDir, kTRUE);
   TString out = outDir + base + TString::Format("_r%d.json", run);

   AtTools::AtCut2D raw(Form("%s_raw_r%d", cut.GetName().c_str(), run), v, cut.GetXAxis(), cut.GetYAxis());
   AtTools::AtParticleID rawPid(raw, pid.GetZ(), pid.GetA());
   if (!rawPid.WriteJSON(out.Data())) {
      std::cout << "\033[1;31mERROR: cannot write " << out << "\033[0m\n";
      return;
   }
   std::cout << "run " << run << ": f=" << f << ", x scaled by 1/sqrt(f)=" << s << " -> " << out << "\n";
}
