/// @file make_sim_points.C
/// @brief Rebuild the truth-matched simulated proton point cache at a chosen fMinPoints.
///
/// gate_draw_sim.C / gate_draw_pp.C read plots/sim_proton_points.root, and the version they build
/// themselves comes from the PERSISTED AtPIDEvent in *_genfitter_pp.root -- which was written at
/// the class default fMinPoints = 30 and cannot be changed without re-running the production. This
/// re-runs AtSpyralPID::Estimate() in memory over the pattern tracks in *_reco.root instead, so the
/// cut is a parameter.
///
/// WHY IT MATTERS FOR A GATE. fMinPoints rejects short tracks, and short means low energy, which
/// means low Brho and high dE/dx -- the far end of the proton locus. Lowering it from 30 to 15 does
/// not move any track already on the plane; it ADDS tracks at exactly the low-Brho end where the
/// gate's lower boundary is drawn and where the theta_cm < 20 deg acceptance is lost. A gate drawn
/// on the 30 plane is drawn on a plane missing the population it most needs to include.
///
/// Truth match is |(180 - polar) - theta_true| < 10 deg: the reconstructed polar is measured
/// against the opposite z sense (the simulation reverses drift z in digitisation).
///
///   root -b -q 'make_sim_points.C(15,"plots/sim_proton_points_mp15.root")'

void make_sim_points(Int_t minPoints = 15, TString outFile = "plots/sim_proton_points_mp15.root",
                     TString simDir = "/mnt/f/a1975_C16_pp_pid",
                     TString tags = "s2001,s2002,s2003,s2004,s2005,s2006", Double_t bField = 2.85,
                     Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) spy.SetMinPoints(minPoints);
   printf("\n  fMinPoints = %d\n", minPoints > 0 ? minPoints : 30);

   std::vector<float> vx, vy;
   long nGen = 0, nMatch = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fs = simDir + "/" + tg + "_sim.root", fr = simDir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tr = Fr ? (TTree *)Fr->Get("cbmsim") : nullptr;
      if (!ts || !tr) { if (Fs) Fs->Close(); if (Fr) Fr->Close(); continue; }
      TClonesArray *mc = nullptr, *pa = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);

      Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i); tr->GetEntry(i);
         double thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp > 0) thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thT < 0) continue;
         ++nGen;
         if (!pa || !pa->GetEntriesFast()) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         double bd = 1e9, bx = 0, by = 0; bool got = false;
         for (auto &track : pe->GetTrackCand()) {
            AtTrack &t2 = const_cast<AtTrack &>(track);
            auto res = spy.Estimate(t2);
            if (!res.valid) continue;
            double d = std::fabs(180.0 - res.polar * TMath::RadToDeg() - thT);
            if (d < bd) { bd = d; bx = res.sqrtdEdx; by = res.brho; got = true; }
         }
         if (got && bd < dThetaMax) { vx.push_back(bx); vy.push_back(by); ++nMatch; }
      }
      Fs->Close(); Fr->Close();
      printf("  %-8s  %ld matched so far\n", tg.Data(), nMatch);
   }
   delete ta;

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   gSystem->mkdir(here + "/plots", kTRUE);
   TFile out(here + "/" + outFile, "RECREATE");
   TTree t("pts", "truth-matched simulated protons");
   float x, y; t.Branch("x", &x); t.Branch("y", &y);
   for (size_t i = 0; i < vx.size(); ++i) { x = vx[i]; y = vy[i]; t.Fill(); }
   t.Write(); out.Close();

   // where the added tracks land: the low-Brho end is the whole reason for doing this
   long lo = 0;
   for (size_t i = 0; i < vy.size(); ++i) if (vy[i] < 0.25) ++lo;
   printf("\n  %zu truth-matched protons of %ld generated (%.1f %%)\n", vx.size(), nGen,
          100.0 * vx.size() / std::max(1L, nGen));
   printf("  below Brho 0.25 T*m: %ld (%.1f %%)\n", lo, 100.0 * lo / std::max(1UL, vx.size()));
   printf("  wrote %s\n\n", outFile.Data());
}
