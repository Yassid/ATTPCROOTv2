/// @file apply_gate_C15p.C
/// @brief Apply a PID gate to the cached plane and write a per-track selection list.
///
///   root -b -q 'apply_gate_C15p.C("gates/proton_C15p.json")'
///   root -b -q 'apply_gate_C15p.C("gates/proton_C15p.json", "/home/yassid/a2091_C15_reco/", "sel_p.root")'
///
/// Reads every <run>_pid.root, tests each track against the polygon, and writes a `sel` tree
/// of (run, event, trackID) for the ones inside. That triple is what makes the selection
/// USABLE: the fitting stage can join on it and fit only the gated tracks, without recomputing
/// AtSpyralPID over the multi-GB recos. A gate you can only draw is half a gate.
///
/// Gate files are spyral_utils Cut2D JSON, written by draw_gate_C15p.C on THIS plane. Applies the
/// same per-run gain table the gate was drawn with, since the cached dE/dx is raw.
///
/// The report breaks the yield down by run so a gate that only holds in part of the run set
/// is visible immediately -- that is the signature of a gain-matching problem, and it would
/// otherwise hide inside a healthy-looking total.

#include "gain_C15p.h"

void apply_gate_C15p(TString gateFile, TString pointsFile = "pid/points_C15p.root", TString outFile = "",
                     Int_t minClusters = 0, Double_t maxVtxR = -1.0, Bool_t perRun = kTRUE,
                     /// MUST match the table the gate was DRAWN on. A gate drawn on a matched
                     /// plane and applied to raw values selects the wrong tracks, silently.
                     /// IC window. MUST MATCH THE ONE THE GATE WAS DRAWN WITH -- a gate drawn on a
                     /// single-beam plane and applied to the full cocktail counts tracks from a beam
                     /// it was never meant to select, and the number looks perfectly reasonable.
                     Double_t icLo = -1, Double_t icHi = -1, Int_t runMin = 138, Int_t runMax = 182,
                     /// The IC window was chosen on the SINGLE-PULSE spectrum and pid/ic_C15p.json
                     /// records singlePulse: true, so applying the window without this condition is
                     /// not the gate that was selected. gate_events_C15p.C requires it too; leaving
                     /// them inconsistent makes the pre-fit and post-fit paths disagree by ~9 %.
                     Bool_t requireSinglePulse = kTRUE)
{
   gSystem->Load("libAtTools.so");

   // gSystem->ExpandPathName has TWO overloads: the char* one returns the expanded path, the
   // TString one expands IN PLACE and returns a Bool_t error flag. Assigning that Bool_t back to
   // the TString blanks it, and the caller then reports "cannot load" against an empty filename.

   gSystem->ExpandPathName(gateFile);
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

   // The points file already carries the gain match and the IC join, with the length checks that
   // join needs. Re-deriving either here would mean two implementations that can disagree.
   if (gSystem->AccessPathName(pointsFile)) {
      std::cout << "\033[1;31mERROR: " << pointsFile << " not found. Build it with "
                   "pid/make_points_C15p.C().\033[0m\n";
      return;
   }
   TFile *fin = TFile::Open(pointsFile);
   TTree *ch = fin ? (TTree *)fin->Get("pts") : nullptr;
   if (!ch) {
      std::cout << "\033[1;31mERROR: no tree 'pts' in " << pointsFile << "\033[0m\n";
      return;
   }

   Int_t run, event, trackID, nclusters, npulse;
   Float_t sqrtdedx, brho, ic, vtxr;
   ch->SetBranchAddress("run", &run);
   ch->SetBranchAddress("event", &event);
   ch->SetBranchAddress("trackID", &trackID);
   ch->SetBranchAddress("nclusters", &nclusters);
   ch->SetBranchAddress("npulse", &npulse);
   ch->SetBranchAddress("sqrtdedx", &sqrtdedx);
   ch->SetBranchAddress("brho", &brho);
   ch->SetBranchAddress("ic", &ic);
   ch->SetBranchAddress("vtxr", &vtxr);

   if (outFile.Length() == 0) {
      TString base = gSystem->BaseName(gateFile);
      base.ReplaceAll(".json", "");
      outFile = TString("pid/sel_") + base + ".root";
   }
   TFile fo(outFile, "RECREATE");
   TTree sel("sel", "tracks inside the PID gate");
   Int_t s_run, s_event, s_track;
   sel.Branch("run", &s_run, "run/I");
   sel.Branch("event", &s_event, "event/I");
   sel.Branch("trackID", &s_track, "trackID/I");

   const bool icOn = (icLo >= 0 && icHi > icLo);
   std::map<int, std::pair<long, long>> perRunCounts;
   Long64_t nAll = 0, nCons = 0, nIn = 0, nNoIC = 0;
   std::set<int> runsSeen;

   for (Long64_t i = 0; i < ch->GetEntries(); ++i) {
      ch->GetEntry(i);
      ++nAll;
      if (run < runMin || run > runMax)
         continue;
      if (minClusters > 0 && nclusters < minClusters)
         continue;
      if (maxVtxR > 0 && vtxr >= maxVtxR)
         continue;
      if (icOn) {
         if (ic < 0) { ++nNoIC; continue; }   // run without usable IC: cannot satisfy a window
         if (ic < icLo || ic > icHi)
            continue;
         if (requireSinglePulse && npulse != 1)
            continue;
      }
      ++nCons;
      runsSeen.insert(run);
      perRunCounts[run].second++;
      if (!cut.IsInside(sqrtdedx, brho))
         continue;
      ++nIn;
      perRunCounts[run].first++;
      s_run = run; s_event = event; s_track = trackID;
      sel.Fill();
   }

   fo.cd();
   sel.Write();
   fo.Close();

   std::cout << "\033[1;33m=== apply_gate_C15p : " << cut.GetName() << " ===\033[0m\n"
             << "  gate     : " << gateFile << "  (" << cut.GetVertices().size() << " vertices)\n"
             << "  runs     : " << runsSeen.size() << "\n"
             << "  tracks   : " << nAll << " total, " << nCons << " considered (valid"
             << (minClusters > 0 ? Form(" && nClusters>=%d", minClusters) : "")
             << (maxVtxR > 0 ? Form(" && vtxR<%g", maxVtxR) : "") << ")\n"
             << "  \033[1;32mIN GATE  : " << nIn << "  = " << (nCons ? 100.0 * nIn / nCons : 0.)
             << "% of considered\033[0m\n"
             << "  selection: " << outFile << "  (run, event, trackID)\n";
   std::cout << "  IC       : " << (icOn ? Form("[%.0f, %.0f]", icLo, icHi) : "OFF (all beam species)")
             << (nNoIC ? Form("  -- %lld tracks dropped for having no usable IC", (long long)nNoIC) : "")
             << "\n";

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
