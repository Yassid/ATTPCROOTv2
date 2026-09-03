/// @file make_points_C15p.C
/// @brief Merge the per-run PID caches into the single `pts` points file the gate GUI reads,
///        applying the per-run gain match and joining the ion chamber onto every track.
///
///   root -b -q 'pid/make_points_C15p.C()'
///   root -b -q 'pid/make_points_C15p.C("/home/yassid/a2091_C15_reco/","pid/points_C15p.root")'
///
/// Output tree `pts`, one entry per VALID track, with the branch names gate_draw_C15p.C expects:
///   sqrtdedx  gain-matched sqrt(dE/dx)          brho   T*m
///   polar     degrees (AtSpyralResult is rad)   ic     ion-chamber max, -1 if unavailable
///   npulse    IC pulse count (pile-up)          run / event / trackID / nclusters
///
/// ★ THE GAIN IS APPLIED HERE, not in the reco, which persists raw dE/dx. So this file is already
/// matched and the GUI must not match it again. The header printed at the end says which table was
/// used; if it says NONE, every gate drawn on the result is drawn on a plane smeared by the ~29 %
/// per-run drift.
///
/// ★ THE IC JOIN IS BY ENTRY INDEX. <run>_ic.root's `entry` is the position in the FRIB tree and
/// the PID cache's `event` is the position in the reco tree; both come from the same HDF5 event
/// ordering, so they correspond. That is only true while the two have the SAME LENGTH, so this
/// checks and refuses to join a run where they differ rather than silently pairing track n with
/// some other event's beam particle -- which would produce a perfectly plausible, entirely wrong
/// IC gate.
///
/// Runs with no IC (empty or absent /frib: 0024, 0025, 0047) get ic = -1 and are counted. A
/// downstream IC window then drops them, which is correct, but it has to be a visible decision
/// rather than a surprise.

#include "../gain_C15p.h"

void make_points_C15p(TString inDir = "/home/yassid/C15p_reco/", TString outFile = "",
                      TString gainTable = "gainmatch_C15p.csv", TString icDir = "/home/yassid/a2091_C15_ic/",
                      Int_t runMin = 138, Int_t runMax = 182)
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (outFile.IsNull())
      outFile = here + "/points_C15p.root";

   TString gt = gainTable;
   if (gt.Length() && !gt.BeginsWith("/") && gSystem->AccessPathName(gt)) {
      gt = here + "/../" + gainTable; // pid/ -> workspace root
      if (gSystem->AccessPathName(gt))
         gt = here + "/" + gainTable;
   }
   auto gain = LoadGainTable_C15p(gt);
   if (gain.empty())
      std::cout << "\033[1;31mWARNING: no gain table -- the points file will be RAW.\033[0m\n";

   TFile fo(outFile, "RECREATE");
   TTree pts("pts", "C15p PID points, gain matched, IC joined");
   Float_t sqrtdedx, brho, polar, ic, dedx, arclen, vtxz, vtxr;
   Int_t run, event, trackID, nclusters, npulse;
   pts.Branch("sqrtdedx", &sqrtdedx, "sqrtdedx/F");
   pts.Branch("brho", &brho, "brho/F");
   pts.Branch("polar", &polar, "polar/F");
   pts.Branch("ic", &ic, "ic/F");
   pts.Branch("npulse", &npulse, "npulse/I");
   pts.Branch("dedx", &dedx, "dedx/F");
   pts.Branch("arclen", &arclen, "arclen/F");
   pts.Branch("vtxz", &vtxz, "vtxz/F");
   pts.Branch("vtxr", &vtxr, "vtxr/F");
   pts.Branch("run", &run, "run/I");
   pts.Branch("event", &event, "event/I");
   pts.Branch("trackID", &trackID, "trackID/I");
   pts.Branch("nclusters", &nclusters, "nclusters/I");

   Long64_t nTot = 0, nNoIC = 0;
   Int_t nRuns = 0, nRunsNoIC = 0, nRunsMismatch = 0, nRunsNoGain = 0;

   for (Int_t r = runMin; r <= runMax; ++r) {
      TString pf = TString::Format("%srun_%04d_pid.root", inDir.Data(), r);
      if (gSystem->AccessPathName(pf))
         continue;
      TFile *fp = TFile::Open(pf);
      TTree *tp = fp ? (TTree *)fp->Get("pid") : nullptr;
      if (!tp) {
         if (fp) fp->Close();
         continue;
      }
      ++nRuns;

      bool missingGain = false;
      const double gf = GainFactor_C15p(gain, r, missingGain);
      if (missingGain)
         ++nRunsNoGain;
      const double sq = std::sqrt(gf);

      // --- IC summary for this run, indexed by entry -------------------------------------
      std::vector<float> icv;
      std::vector<int> npv;
      TString icf = TString::Format("%srun_%04d_ic.root", icDir.Data(), r);
      if (!gSystem->AccessPathName(icf)) {
         TFile *fi = TFile::Open(icf);
         TTree *ti = fi ? (TTree *)fi->Get("ic") : nullptr;
         if (ti) {
            Int_t e_, np_;
            Float_t im_;
            ti->SetBranchAddress("entry", &e_);
            ti->SetBranchAddress("icmax", &im_);
            ti->SetBranchAddress("npulse", &np_);
            const Long64_t ni = ti->GetEntries();
            icv.assign(ni, -1.f);
            npv.assign(ni, 0);
            for (Long64_t i = 0; i < ni; ++i) {
               ti->GetEntry(i);
               if (e_ >= 0 && e_ < (Int_t)icv.size()) {
                  icv[e_] = im_;
                  npv[e_] = np_;
               }
            }
         }
         if (fi) fi->Close();
      }
      if (icv.empty())
         ++nRunsNoIC;

      Int_t b_run, b_event, b_track, b_valid, b_nc, b_dir, b_np;
      Double_t b_sq, b_brho, b_pol, b_dedx, b_arc, b_vz, b_vr;
      tp->SetBranchAddress("run", &b_run);
      tp->SetBranchAddress("event", &b_event);
      tp->SetBranchAddress("trackID", &b_track);
      tp->SetBranchAddress("valid", &b_valid);
      tp->SetBranchAddress("nClusters", &b_nc);
      tp->SetBranchAddress("sqrtdEdx", &b_sq);
      tp->SetBranchAddress("brho", &b_brho);
      tp->SetBranchAddress("polar", &b_pol);
      tp->SetBranchAddress("dEdx", &b_dedx);
      tp->SetBranchAddress("arclength", &b_arc);
      tp->SetBranchAddress("vtxZ", &b_vz);
      tp->SetBranchAddress("vtxR", &b_vr);

      // ★ THE JOIN IS POSITIONAL, SO THE TWO TREES MUST HAVE THE SAME LENGTH -- and "the index
      // fits" is NOT the test. Measured on this dataset, FRIB and GET agree exactly on most runs
      // but not all: run_0023 has 360 FRIB events against 18,862 GET (caught by any range test),
      // and run_0022 has 58,840 against 39,978 -- MORE. A too-long IC file passes an in-range
      // check and pairs every track with some other event's beam particle, producing a perfectly
      // plausible and entirely wrong IC gate. So require the counts to MATCH.
      //
      // The reference count is the reco tree's entry count, read from the file header (cheap, no
      // event data). If the reco is unavailable, fall back to the highest event index seen, which
      // undercounts only by however many trailing events carried no valid track.
      Long64_t nReco = -1;
      TString rf = TString::Format("%srun_%04d_reco.root", inDir.Data(), r);
      if (!gSystem->AccessPathName(rf)) {
         TFile *fr2 = TFile::Open(rf);
         if (fr2 && !fr2->IsZombie()) {
            if (auto *tr = (TTree *)fr2->Get("cbmsim"))
               nReco = tr->GetEntries();
            fr2->Close();
         }
      }
      Long64_t maxEvent = -1;
      for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
         tp->GetEntry(i);
         if (b_event > maxEvent)
            maxEvent = b_event;
      }
      const Long64_t nRef = nReco > 0 ? nReco : maxEvent + 1;
      const bool exact = nReco > 0;
      bool useIC = !icv.empty();
      if (useIC) {
         const double rel = nRef > 0 ? std::abs((double)icv.size() - nRef) / nRef : 1.0;
         // Exact equality when the reco count is known; 2 % otherwise, which covers only the
         // trailing-events-without-tracks slack in the fallback.
         // ±1 event is a boundary artifact -- one stream stopping an event short of the other --
         // and indices 0..min-1 still line up, so it is safe to join over the overlap. ANYTHING
         // LARGER IS REFUSED, even 61 events on 25,000 (0.24 %): a surplus that is not at the end
         // but interleaved shifts every later event, and there is no way to tell the two apart
         // without event IDs, which the IC summary does not carry. Being wrong here selects the
         // wrong beam silently, so the tolerance stays at what can be explained.
         const Long64_t diff = std::abs((Long64_t)icv.size() - nRef);
         const bool ok = exact ? (diff <= 1) : (rel <= 0.02);
         if (!ok) {
            std::cout << "\033[1;31m  run " << r << ": IC has " << icv.size() << " entries vs " << nRef
                      << (exact ? " reco events" : " (est.) events")
                      << " -- REFUSING to join, ic set to -1.\033[0m\n";
            useIC = false;
            ++nRunsMismatch;
         }
      }

      for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
         tp->GetEntry(i);
         if (b_valid != 1)
            continue;
         run = b_run;
         event = b_event;
         trackID = b_track;
         nclusters = b_nc;
         sqrtdedx = b_sq * sq;
         dedx = b_dedx * gf;
         brho = b_brho;
         polar = b_pol * TMath::RadToDeg();
         arclen = b_arc;
         vtxz = b_vz;
         vtxr = b_vr;
         if (useIC && b_event >= 0 && b_event < (Long64_t)icv.size()) {
            ic = icv[b_event];
            npulse = npv[b_event];
         } else {
            ic = -1;
            npulse = 0;
            ++nNoIC;
         }
         pts.Fill();
         ++nTot;
      }
      fp->Close();
   }

   fo.cd();
   pts.Write();
   fo.Close();

   std::cout << "\033[1;33m=== make_points_C15p ===\033[0m\n"
             << "  runs        : " << nRuns << "  (" << runMin << "-" << runMax << ")\n"
             << "  tracks      : " << nTot << " valid\n"
             << "  gain        : " << (gain.empty() ? "NONE (RAW)" : gt.Data()) << "\n"
             << "  runs w/o IC : " << nRunsNoIC << (nRunsMismatch ? Form(", %d refused on length mismatch", nRunsMismatch) : "")
             << "  -> " << nNoIC << " tracks carry ic = -1\n";
   if (nRunsNoGain)
      std::cout << "\033[1;31m  " << nRunsNoGain << " run(s) had no gain entry and went in UNMATCHED\033[0m\n";
   std::cout << "  \033[1;32mwrote\033[0m " << outFile << "\n";
}
