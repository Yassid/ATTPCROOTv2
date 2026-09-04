/// @file make_points_C15d.C
/// @brief Merge the per-run PID caches into the single `pts` points file the gate GUI reads,
///        applying the per-run gain match and joining the ion chamber onto every track.
///
///   root -b -q 'pid/make_points_C15d.C()'
///   root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root")'
///
/// Output tree `pts`, one entry per VALID track, with the branch names gate_draw_C15d.C expects:
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

#include "../gain_C15d.h"

void make_points_C15d(TString inDir = "/home/yassid/C15d_reco/", TString outFile = "",
                      TString gainTable = "gainmatch_C15d.csv", TString icDir = "/home/yassid/C15d_ic/",
                      Int_t runMin = 13, Int_t runMax = 133)
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (outFile.IsNull())
      outFile = here + "/points_C15d.root";

   TString gt = gainTable;
   if (gt.Length() && !gt.BeginsWith("/") && gSystem->AccessPathName(gt)) {
      gt = here + "/../" + gainTable; // pid/ -> workspace root
      if (gSystem->AccessPathName(gt))
         gt = here + "/" + gainTable;
   }
   auto gain = LoadGainTable_C15d(gt);
   if (gain.empty())
      std::cout << "\033[1;31mWARNING: no gain table -- the points file will be RAW.\033[0m\n";

   TFile fo(outFile, "RECREATE");
   TTree pts("pts", "C15d PID points, gain matched, IC joined");
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
   Int_t nRuns = 0, nRunsNoIC = 0, nRunsMismatch = 0, nRunsNoGain = 0, nRunsTruncated = 0;

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
      const double gf = GainFactor_C15d(gain, r, missingGain);
      if (missingGain)
         ++nRunsNoGain;
      const double sq = std::sqrt(gf);

      // --- IC summary for this run, indexed by entry -------------------------------------
      std::vector<float> icv;
      std::vector<int> npv;
      bool hasEvtId = false, icGapless = false;
      Long64_t nIcEntries = 0;
      TString icf = TString::Format("%srun_%04d_ic.root", icDir.Data(), r);
      if (!gSystem->AccessPathName(icf)) {
         TFile *fi = TFile::Open(icf);
         TTree *ti = fi ? (TTree *)fi->Get("ic") : nullptr;
         if (ti) {
            Int_t e_, np_, ev_ = -1;
            Float_t im_;
            ti->SetBranchAddress("entry", &e_);
            ti->SetBranchAddress("icmax", &im_);
            ti->SetBranchAddress("npulse", &np_);
            // ★ evtid is the TRUE event number, taken from the HDF5 dataset name (evt<N>_1903).
            // It is written only by the corrected icsum_C15d.C. OLDER IC SUMMARIES DO NOT HAVE IT
            // and must keep the previous behaviour exactly -- other experiments (16C, 12Be) join
            // cleanly on position and nothing here may change for them.
            hasEvtId = (ti->GetBranch("evtid") != nullptr);
            if (hasEvtId)
               ti->SetBranchAddress("evtid", &ev_);
            const Long64_t ni = ti->GetEntries();
            // Size the arrays by the LARGEST event number seen, not by the entry count: with
            // evtid the IC is indexed by event number, so a run whose IC stops early still places
            // its events at the right indices.
            Long64_t maxIdx = ni;
            if (hasEvtId) {
               for (Long64_t i = 0; i < ni; ++i) { ti->GetEntry(i); if (ev_ >= 0 && ev_ + 1 > maxIdx) maxIdx = ev_ + 1; }
               // ★ SKIP EMPTY EVENTS. An event with no FRIB payload is written with evtid = -1
               // (and npulse 0, icmax -1). Counting that as a gap is what refused the join on
               // every run: run_0027 has entry == evtid for all 54731 real events and exactly ONE
               // trailing empty one, and that single entry was enough to discard the run.
               icGapless = true;
               Int_t prev = -1;
               for (Long64_t i = 0; i < ni; ++i) {
                  ti->GetEntry(i);
                  if (ev_ < 0) continue;                       // empty event, not a gap
                  if (prev >= 0 && ev_ != prev + 1) { icGapless = false; break; }
                  prev = ev_;
               }
            }
            icv.assign(maxIdx, -1.f);
            npv.assign(maxIdx, 0);
            for (Long64_t i = 0; i < ni; ++i) {
               ti->GetEntry(i);
               const Int_t idx = hasEvtId ? ev_ : e_;
               if (idx >= 0 && idx < (Int_t)icv.size()) {
                  icv[idx] = im_;
                  npv[idx] = np_;
               }
            }
            nIcEntries = ni;
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
         bool ok = exact ? (diff <= 1) : (rel <= 0.02);
         // ★ NEW PATH, and it only ever ACCEPTS what the old code would have refused -- it can
         // never reject something the old code accepted, so 16C/12Be are bit-for-bit unaffected.
         // With evtid present AND the sequence gapless, the IC is indexed by TRUE EVENT NUMBER,
         // so a size difference means the IC simply STOPPED EARLY: events 0..n_ic-1 match exactly
         // and the rest legitimately have no IC. Verified on a2091 run_0027: evtid 0..45863 with
         // ZERO gaps against 55006 pad events. Without evtid, or with gaps, fall through to the
         // old tolerance -- an interleaved surplus shifts every later event and is unrecoverable.
         if (!ok && hasEvtId && icGapless) {
            // Report the direction correctly. EITHER side can be shorter: the IC can stop early
            // (reco events past its end get ic = -1) or the RECO can be short, which means the
            // event range was under-derived and that run should be re-reconstructed rather than
            // quietly analysed at partial statistics. Printing a negative count hid the second
            // case entirely.
            const Long64_t d = nRef - (Long64_t)nIcEntries;
            if (d >= 0)
               std::cout << "  run " << r << ": IC has " << nIcEntries << " events vs " << nRef
                         << " reco -- gapless evtid, joining over the overlap (" << d
                         << " reco events beyond the IC get ic = -1)\n";
            else
               // NOT an error. /get (pad plane) and /frib (auxiliary) are separate DAQs and can
               // record different numbers of events in the same run -- run_0118 has 31446 GET
               // against 64844 FRIB, run_0098 has 1 against 150893. The reconstruction is complete
               // for what /get holds; the surplus FRIB events simply have no pad data. Reported so
               // the asymmetry is visible, not as a fault.
               std::cout << "  run " << r << ": IC has " << nIcEntries << " events, reco has " << nRef
                         << " (" << (100.0 * nRef / nIcEntries) << " % -- fewer GET than FRIB events, normal);"
                         << " joining the overlap\n";
            ok = true;
            ++nRunsTruncated;
         }
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

   std::cout << "\033[1;33m=== make_points_C15d ===\033[0m\n"
             << "  runs        : " << nRuns << "  (" << runMin << "-" << runMax << ")\n"
             << "  tracks      : " << nTot << " valid\n"
             << "  gain        : " << (gain.empty() ? "NONE (RAW)" : gt.Data()) << "\n"
             << "  runs w/o IC : " << nRunsNoIC << (nRunsMismatch ? Form(", %d refused on length mismatch", nRunsMismatch) : "")
             << (nRunsTruncated ? Form(", %d joined over the overlap via evtid", nRunsTruncated) : "")
             << "  -> " << nNoIC << " tracks carry ic = -1\n";
   if (nRunsNoGain)
      std::cout << "\033[1;31m  " << nRunsNoGain << " run(s) had no gain entry and went in UNMATCHED\033[0m\n";
   std::cout << "  \033[1;32mwrote\033[0m " << outFile << "\n";
}
