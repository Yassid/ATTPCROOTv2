/// @file gate_events_C15d.C
/// @brief Build a small FIT INPUT: only IC-gated events, and inside them only PID-gated tracks.
///
///   root -b -q 'pid/gate_events_C15d.C("run_0026","pid/proton_C15d.json")'
///
/// Same construction as a1954's gate_events_C14.C -- clone the reco tree, keep the events that
/// pass the beam window, replace each event's track list with the gated subset, and write the
/// FairRoot metadata (BranchList / TimeBasedBranchList / FileHeader / cbmout) so the fitters can
/// read the result. Two things differ, both because of how this workspace stores its PID.
///
/// ★ GATE BEFORE FITTING. Fitting is the expensive stage, and the gate exists before any fit is
/// needed. On the (d,d') pass, fitting ungated meant ~1.9M fits to keep ~133k deuterons.
///
/// ★ THE GATE IS TESTED ON GAIN-MATCHED dE/dx. Every gate here is drawn on the matched plane while
/// AtSpyralPID returns raw values, so the per-run factor is applied before the polygon test:
///
///     sqrtdEdx_matched = sqrtdEdx_raw * sqrt(f(run))
///
/// Handing the matched gate to a raw quantity selects a different set of tracks. That is not
/// hypothetical: passing a per-run RAW-rescaled gate to AtGenfitter::SetPIDGate instead selected
/// 4,217 tracks on run_0026 where this selection gives 2,606 -- 62 % more, silently. The fitter's
/// internal AtSpyralPID is not the same instance that built the plane, and does not have to agree
/// with it, which is why the selection is made here against the SAME numbers the plane was drawn
/// on rather than delegated to the fitter.
///
/// It uses the PERSISTED AtPIDEvent when present (that is literally the plane), matching tracks by
/// trackID, and only falls back to recomputing AtSpyralPID if the branch is missing.

#include "../gain_C15d.h"

void gate_events_C15d(TString run, TString gateFile = "pid/proton_C15d.json",
                      TString inDir = "/home/yassid/C15d_reco/",
                      TString outDir = "/home/yassid/C15d_fit/in/", TString icDir = "/home/yassid/C15d_ic/",
                      Double_t icLo = 931, Double_t icHi = 1413, TString gainTable = "gainmatch_C15d.csv",
                      Double_t bField = 2.85, Double_t thMinDeg = -1, Bool_t requireSinglePulse = kTRUE)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const Int_t runNo = AtGainMatchTask::RunNumberFromName(run);

   TString gt = gainTable;
   if (gt.Length() && !gt.BeginsWith("/") && gSystem->AccessPathName(gt)) {
      gt = here + "/../" + gainTable;
      if (gSystem->AccessPathName(gt))
         gt = here + "/" + gainTable;
   }
   auto gain = LoadGainTable_C15d(gt, false);
   bool missingGain = false;
   const double gf = GainFactor_C15d(gain, runNo, missingGain);
   if (missingGain) {
      printf("\033[1;31m%s: no gain entry -- refusing to gate, the polygon would be tested against "
             "unmatched dE/dx.\033[0m\n", run.Data());
      return;
   }
   const double sq = std::sqrt(gf);

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   const auto &cut = pid.GetCut();
   if (!cut.IsValid()) {
      printf("\033[1;31mERR: cannot load gate %s\033[0m\n", gateFile.Data());
      return;
   }

   TString rf = inDir + run + "_reco.root";
   if (gSystem->AccessPathName(rf)) {
      printf("SKIP %s (no reco)\n", run.Data());
      return;
   }

   // ---- IC keep-flags, with the same length check make_points_C15d.C uses --------------------
   TFile *fin = TFile::Open(rf);
   TTree *t = fin ? (TTree *)fin->Get("cbmsim") : nullptr;
   if (!t) {
      printf("SKIP %s (no tree)\n", run.Data());
      if (fin) fin->Close();
      return;
   }
   const Long64_t nReco = t->GetEntries();

   std::vector<char> icKeep(nReco, 0);
   bool icOn = (icLo >= 0 && icHi > icLo);
   Long64_t nIC = 0;
   if (icOn) {
      TString icf = icDir + run + "_ic.root";
      TFile *fi = gSystem->AccessPathName(icf) ? nullptr : TFile::Open(icf);
      TTree *ti = fi ? (TTree *)fi->Get("ic") : nullptr;
      if (!ti) {
         printf("\033[1;33m%s: no IC summary -- refusing to write an ungated file under a gated "
                "name.\033[0m\n", run.Data());
         if (fi) fi->Close();
         fin->Close();
         return;
      }
      // FRIB and GET agree on most runs but not all (run_0022 has 47 % MORE FRIB events, run_0023
      // 98 % fewer). A positional join across a mismatch pairs tracks with another event's beam.
      if (std::llabs(ti->GetEntries() - nReco) > 1) {
         printf("\033[1;31m%s: IC has %lld entries vs %lld reco events -- REFUSING to gate.\033[0m\n",
                run.Data(), (long long)ti->GetEntries(), (long long)nReco);
         fi->Close();
         fin->Close();
         return;
      }
      Int_t e_, np_;
      Float_t im_;
      ti->SetBranchAddress("entry", &e_);
      ti->SetBranchAddress("icmax", &im_);
      ti->SetBranchAddress("npulse", &np_);
      for (Long64_t i = 0; i < ti->GetEntries(); ++i) {
         ti->GetEntry(i);
         if (e_ < 0 || e_ >= (Int_t)icKeep.size())
            continue;
         if (im_ >= icLo && im_ <= icHi && (!requireSinglePulse || np_ == 1)) {
            icKeep[e_] = 1;
            ++nIC;
         }
      }
      fi->Close();
   } else {
      std::fill(icKeep.begin(), icKeep.end(), 1);
      nIC = nReco;
   }

   // ---- filter ------------------------------------------------------------------------------
   TClonesArray *pe = nullptr, *pidArr = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);
   const bool hasPID = t->GetBranch("AtPIDEvent") != nullptr;
   if (hasPID)
      t->SetBranchAddress("AtPIDEvent", &pidArr);
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);
   if (!hasPID)
      printf("\033[1;33m%s: no AtPIDEvent branch -- recomputing AtSpyralPID (slower, and it must "
             "match the settings the plane was built with)\033[0m\n", run.Data());

   gSystem->mkdir(outDir.Data(), kTRUE);
   TString of = outDir + run + "_reco.root";
   TFile *fout = new TFile(of, "RECREATE", "", 1);
   TTree *nt = t->CloneTree(0);

   Long64_t nEvt = 0, nTrk = 0, nSeen = 0;
   for (Long64_t i = 0; i < nReco; ++i) {
      if (!icKeep[i])
         continue;
      t->GetEntry(i);
      if (!pe || pe->GetEntries() == 0)
         continue;
      auto *p = (AtPatternEvent *)pe->At(0);
      if (!p)
         continue;

      // trackID -> matched (sqrtdEdx, brho) from the persisted plane
      std::map<int, std::pair<double, double>> plane;
      if (hasPID && pidArr && pidArr->GetEntriesFast() > 0) {
         auto *pev = (AtPIDEvent *)pidArr->At(0);
         if (pev)
            for (const auto &r : pev->GetSpyral())
               if (r.valid)
                  plane[r.trackID] = {r.sqrtdEdx * sq, r.brho};
      }

      std::vector<AtTrack> keep;
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         ++nSeen;
         double x, y, polDeg;
         auto it = plane.find(tr.GetTrackID());
         if (it != plane.end()) {
            x = it->second.first;
            y = it->second.second;
            polDeg = -1; // polar not needed unless thMin is set; recomputed below if it is
         } else if (!hasPID) {
            auto r = spy.Estimate(tr);
            if (!r.valid)
               continue;
            x = r.sqrtdEdx * sq;
            y = r.brho;
            polDeg = r.polar * TMath::RadToDeg();
         } else {
            continue; // track had no valid PID entry: it is not on the plane the gate was drawn on
         }
         if (thMinDeg > 0) {
            if (polDeg < 0) {
               auto r = spy.Estimate(tr);
               polDeg = r.valid ? r.polar * TMath::RadToDeg() : -1;
            }
            if (!(polDeg > thMinDeg))
               continue;
         }
         if (!cut.IsInside(x, y))
            continue;
         keep.push_back(tr);
      }
      if (keep.empty())
         continue;
      nTrk += keep.size();
      p->SetTrackCand(std::move(keep));
      nt->Fill();
      ++nEvt;
   }
   nt->Write();

   // FairRoot metadata, or the fitters cannot open the result
   TList bl;
   for (auto b : *nt->GetListOfBranches())
      bl.Add(new TObjString(b->GetName()));
   fout->cd();
   bl.Write("BranchList", TObject::kSingleKey);
   TList etb;
   etb.Write("TimeBasedBranchList", TObject::kSingleKey);
   {
      TFile *fr = TFile::Open(rf, "READ");
      if (fr && !fr->IsZombie()) {
         if (auto *fh = fr->Get("FileHeader")) {
            fout->cd();
            fh->Write("FileHeader");
         }
         if (auto *cb = fr->Get("cbmout")) {
            fout->cd();
            cb->Write("cbmout", TObject::kSingleKey);
         }
         fr->Close();
      }
   }
   fout->Write("", TObject::kOverwrite);
   printf("\033[1;32m%s\033[0m: %lld/%lld events pass IC, %lld tracks seen -> %lld events, %lld tracks "
          "kept (f=%.4f) -> %s\n",
          run.Data(), nIC, nReco, nSeen, nEvt, nTrk, gf, of.Data());
   fout->Close();
   fin->Close();
}
