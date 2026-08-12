/// @file make_data_points.C
/// @brief Rebuild the DATA PID plane at a chosen fMinPoints, for gate drawing.
///
/// gate_draw_pp.C builds its data plane from the persisted AtPIDEvent in *_genfitter_pphand.root,
/// which AtPIDTask wrote at the class default fMinPoints = 30. That cut cannot be changed without
/// re-running the production, so this re-runs AtSpyralPID::Estimate() in memory over the pattern
/// tracks in *_reco.root instead, exactly as AtPIDTask::Exec drives it.
///
/// This is safe to do precisely because the PID does not depend on the fit -- verified on s2001,
/// three paths (in-memory / AtPIDTask standalone / AtPIDTask inside a fit job), 5101 tracks, zero
/// difference in Brho to machine precision. So a plane rebuilt here is the same plane the
/// production gated on, at a different fMinPoints.
///
/// WHY 15 RATHER THAN 30. fMinPoints rejects short tracks, and short means low energy, low Brho,
/// high dE/dx -- the far end of the locus, which is where the gate's lower boundary is drawn and
/// where the theta_cm < 20 deg acceptance is currently lost. On the simulation, going 30 -> 15 adds
/// 2097 protons of which 356 sit below Brho 0.25. Drawing the boundary on a 30 plane means drawing
/// it on a plane missing the tracks it is meant to include.
///
/// NO IC GATE AND NO PID GATE, matching what gate_draw_pp.C shows: every pattern track that
/// produces a valid Spyral result goes on the plane. Filtering here would hide the contamination
/// the gate exists to exclude.
///
///   root -b -q 'pid/make_data_points.C(15,"pid/plots/data_points_mp15.root",106,145,8000)'

void make_data_points(Int_t minPoints = 15, TString outFile = "pid/plots/data_points_mp15.root", Int_t runLo = 106,
                      Int_t runHi = 145, Long64_t maxEvtPerRun = 8000, TString inDir = "/mnt/f/a1975/reco/",
                      Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) spy.SetMinPoints(minPoints);
   printf("\n  fMinPoints = %d, runs %d-%d, <= %lld events each\n\n", minPoints > 0 ? minPoints : 30, runLo, runHi,
          maxEvtPerRun);

   std::vector<float> vx, vy;
   long nTrk = 0, nRun = 0;

   for (int r = runLo; r <= runHi; ++r) {
      TString fn = TString::Format("%srun_%04d_reco.root", inDir.Data(), r);
      if (gSystem->AccessPathName(fn)) continue;
      TFile *f = TFile::Open(fn);
      TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t) { if (f) f->Close(); continue; }
      // the reco file also carries the corrected-event branch; only the pattern is needed
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      TClonesArray *pa = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pa);
      Long64_t N = maxEvtPerRun > 0 ? std::min(t->GetEntries(), maxEvtPerRun) : t->GetEntries();
      long before = (long)vx.size();
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (!pa || !pa->GetEntriesFast()) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         for (auto &track : pe->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(track);
            ++nTrk;
            auto res = spy.Estimate(tr);
            if (!res.valid) continue;
            vx.push_back(res.sqrtdEdx);
            vy.push_back(res.brho);
         }
      }
      f->Close();
      ++nRun;
      printf("  run_%04d  +%ld  (total %zu)\n", r, (long)vx.size() - before, vx.size());
   }

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString full = here + "/../" + outFile;
   gSystem->mkdir(gSystem->DirName(full), kTRUE);
   TFile out(full, "RECREATE");
   TTree t("pts", "data PID points");
   float x, y;
   t.Branch("x", &x);
   t.Branch("y", &y);
   for (size_t i = 0; i < vx.size(); ++i) { x = vx[i]; y = vy[i]; t.Fill(); }
   t.Write();
   out.Close();

   long lo = 0;
   for (size_t i = 0; i < vy.size(); ++i) if (vy[i] < 0.25) ++lo;
   printf("\n  %ld runs, %ld pattern tracks, %zu valid PID points (%.1f %%)\n", nRun, nTrk, vx.size(),
          100.0 * vx.size() / std::max(1L, nTrk));
   printf("  below Brho 0.25 T*m: %ld (%.1f %%)\n", lo, 100.0 * lo / std::max(1UL, vx.size()));
   printf("  wrote %s\n\n", outFile.Data());
}
