/// @file apply_gate_C15d.C
/// @brief Apply a PID gate to the cached plane and write a per-track selection list.
///
///   root -b -q 'apply_gate_C15d.C("gates/proton_C15d.json")'
///   root -b -q 'apply_gate_C15d.C("gates/proton_C15d.json", "/home/yassid/C15d_reco/", "sel_p.root")'
///
/// Reads every <run>_pid.root, tests each track against the polygon, and writes a `sel` tree
/// of (run, event, trackID) for the ones inside. That triple is what makes the selection
/// USABLE: the fitting stage can join on it and fit only the gated tracks, without recomputing
/// AtSpyralPID over the multi-GB recos. A gate you can only draw is half a gate.
///
/// Gate files are spyral_utils Cut2D JSON, written by draw_gate_C15d.C on THIS plane. Applies the
/// same per-run gain table the gate was drawn with, since the cached dE/dx is raw.
///
/// The report breaks the yield down by run so a gate that only holds in part of the run set
/// is visible immediately -- that is the signature of a gain-matching problem, and it would
/// otherwise hide inside a healthy-looking total.

#include "gain_C15d.h"

void apply_gate_C15d(TString gateFile, TString inDir = "/home/yassid/C15d_reco/", TString outFile = "",
                     Int_t minClusters = 0, Double_t maxVtxR = -1.0, Bool_t perRun = kTRUE,
                     /// MUST match the table the gate was DRAWN on. A gate drawn on a matched
                     /// plane and applied to raw values selects the wrong tracks, silently.
                     TString gainTable = "gainmatch_C15d.csv", Int_t runMin = 17, Int_t runMax = 103)
{
   gSystem->Load("libAtTools.so");

   gateFile = gSystem->ExpandPathName(gateFile);
   auto cut = AtTools::AtCut2D::LoadJSON(gateFile.Data());
   if (!cut.IsValid()) {
      std::cout << "\033[1;31mERROR: cannot load a valid gate from " << gateFile << "\033[0m\n";
      return;
   }
   if (cut.GetXAxis() != "sqrt_dEdx" && cut.GetXAxis() != "sqrtdEdx")
      std::cout << "\033[1;33mWARNING: gate x-axis is '" << cut.GetXAxis()
                << "', not sqrt_dEdx -- applying it to sqrtdEdx anyway. Check this is intended.\033[0m\n";
   if (cut.GetYAxis() != "brho")
      std::cout << "\033[1;33mWARNING: gate y-axis is '" << cut.GetYAxis() << "', not brho.\033[0m\n";

   TChain ch("pid");
   TSystemDirectory dir(inDir, inDir);
   TList *files = dir.GetListOfFiles();
   if (files == nullptr) {
      std::cout << "\033[1;31mERROR: cannot list " << inDir << "\033[0m\n";
      return;
   }
   Int_t nRuns = 0;
   TIter next(files);
   while (auto *o = dynamic_cast<TSystemFile *>(next())) {
      TString name = o->GetName();
      if (!o->IsDirectory() && name.EndsWith("_pid.root")) {
         TString d = name;
         d.ReplaceAll("run_", "");
         d.ReplaceAll("_pid.root", "");
         const Int_t rn = d.Atoi();
         if (rn < runMin || rn > runMax)
            continue;
         ch.Add(inDir + name);
         ++nRuns;
      }
   }
   if (nRuns == 0) {
      std::cout << "\033[1;31mERROR: no *_pid.root in " << inDir << "\033[0m\n";
      return;
   }

   Int_t run, event, trackID, valid, nClusters;
   Double_t sqrtdEdx, brho, vtxR;
   ch.SetBranchAddress("run", &run);
   ch.SetBranchAddress("event", &event);
   ch.SetBranchAddress("trackID", &trackID);
   ch.SetBranchAddress("valid", &valid);
   ch.SetBranchAddress("nClusters", &nClusters);
   ch.SetBranchAddress("sqrtdEdx", &sqrtdEdx);
   ch.SetBranchAddress("brho", &brho);
   ch.SetBranchAddress("vtxR", &vtxR);

   if (outFile.Length() == 0) {
      TString base = gSystem->BaseName(gateFile);
      base.ReplaceAll(".json", "");
      outFile = inDir + "sel_" + base + ".root";
   }
   TFile fo(outFile, "RECREATE");
   TTree sel("sel", "tracks inside the PID gate");
   Int_t s_run, s_event, s_track;
   sel.Branch("run", &s_run, "run/I");
   sel.Branch("event", &s_event, "event/I");
   sel.Branch("trackID", &s_track, "trackID/I");

   auto gainMap = LoadGainTable_C15d(gainTable);
   Long64_t nMissingGain = 0;

   std::map<int, std::pair<long, long>> perRunCounts; // run -> (in, considered)
   Long64_t nAll = 0, nCons = 0, nIn = 0;

   const Long64_t nEnt = ch.GetEntries();
   for (Long64_t i = 0; i < nEnt; ++i) {
      ch.GetEntry(i);
      ++nAll;
      if (valid != 1)
         continue;
      if (minClusters > 0 && nClusters < minClusters)
         continue;
      if (maxVtxR > 0 && vtxR >= maxVtxR)
         continue;
      ++nCons;
      perRunCounts[run].second++;
      bool missing = false;
      const double gf = GainFactor_C15d(gainMap, run, missing);
      if (missing)
         ++nMissingGain;
      const double sqrtdEdxMatched = sqrtdEdx * std::sqrt(gf);
      if (!cut.IsInside(sqrtdEdxMatched, brho))
         continue;
      ++nIn;
      perRunCounts[run].first++;
      s_run = run;
      s_event = event;
      s_track = trackID;
      sel.Fill();
   }

   fo.cd();
   sel.Write();
   fo.Close();

   std::cout << "\033[1;33m=== apply_gate_C15d : " << cut.GetName() << " ===\033[0m\n"
             << "  gate     : " << gateFile << "  (" << cut.GetVertices().size() << " vertices)\n"
             << "  runs     : " << nRuns << "\n"
             << "  tracks   : " << nAll << " total, " << nCons << " considered (valid"
             << (minClusters > 0 ? Form(" && nClusters>=%d", minClusters) : "")
             << (maxVtxR > 0 ? Form(" && vtxR<%g", maxVtxR) : "") << ")\n"
             << "  \033[1;32mIN GATE  : " << nIn << "  = " << (nCons ? 100.0 * nIn / nCons : 0.)
             << "% of considered\033[0m\n"
             << "  selection: " << outFile << "  (run, event, trackID)\n";
   if (nMissingGain > 0)
      std::cout << "\033[1;31m  WARNING: " << nMissingGain
                << " considered tracks had no gain-table entry and were gated UNMATCHED.\033[0m\n";

   if (perRun && !perRunCounts.empty()) {
      std::cout << "\n  per-run yield (a gate that only holds over part of the run set is a "
                   "gain-matching problem, not a physics one):\n";
      double lo = 1e30, hi = -1e30;
      int loRun = -1, hiRun = -1;
      for (const auto &kv : perRunCounts) {
         if (kv.second.second < 50)
            continue; // too few tracks to mean anything
         const double f = 100.0 * kv.second.first / kv.second.second;
         if (f < lo) { lo = f; loRun = kv.first; }
         if (f > hi) { hi = f; hiRun = kv.first; }
      }
      if (loRun > 0)
         std::cout << "    in-gate fraction spans " << lo << "% (run " << loRun << ") to " << hi << "% (run "
                   << hiRun << ") over runs with >=50 tracks\n";
   }
}
