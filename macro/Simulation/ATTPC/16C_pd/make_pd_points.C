/// @file make_pd_points.C
/// @brief Truth-matched simulated DEUTERON points, for drawing a (p,d) PID gate on the sim plane.
///
/// Reads the PERSISTED AtPIDEvent out of <state>_genfitter.root rather than re-running Estimate()
/// over the reco files. That is both faster (100 MB per state against 1.9 GB) and the RIGHT choice
/// here: the whole (p,d) chain -- data production included -- runs at AtSpyralPID's default
/// fMinPoints = 30, which is exactly what the persisted branch was written at. Re-estimating would
/// only be needed to change that threshold, and changing it is not on the table for (p,d).
///
/// The (p,d) simulation was fitted UNGATED (the gate argument is "" in every <state>_fit.log), so
/// AtPIDEvent carries every pattern track, not just fitted ones -- the plane is complete.
///
/// Truth match: |(180 - polar) - theta_true| < dThetaMax, on the deuteron (PDG 1000010020). The
/// reconstructed polar comes out as 180 - theta_true because the simulation reverses drift z in
/// digitisation; matching polar directly gives exactly zero.
///
/// WHY DRAW A SIM GATE AT ALL. deuteron_tight.json, drawn on the DATA plane, keeps only 79.8 % of
/// truth-matched simulated deuterons, and the loss is state-dependent -- 75.7 % (gs), 77.2 (ex1),
/// 80.5 (ex2), 85.9 (ex3). The rejected tracks form a coherent band along the polygon's upper-left
/// edge, i.e. the same locus sitting slightly high in Brho, not background. Using the data gate on
/// the simulation would fold that 10-point state-dependent efficiency into the acceptance and
/// penalise the ground state most.
///
///   root -b -q 'make_pd_points.C()'

void make_pd_points(TString outFile = "diagnostics/sim_deuteron_points.root",
                    TString simDir = "/mnt/f/a1975_C16_pd_sim",
                    TString states = "gs_s1001,ex1_s1001,ex2_s1001,ex3_s1001", Int_t pdg = 1000010020,
                    Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");

   std::vector<float> vx, vy;
   long nGen = 0;

   TObjArray *ta = states.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fs = simDir + "/" + tg + "_sim.root", ff = simDir + "/" + tg + "_genfitter.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(ff)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *Fs = TFile::Open(fs), *Ff = TFile::Open(ff);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tf = Ff ? (TTree *)Ff->Get("cbmsim") : nullptr;
      if (!ts || !tf) { if (Fs) Fs->Close(); if (Ff) Ff->Close(); continue; }
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchStatus("*", 0);
      tf->SetBranchStatus("AtPIDEvent*", 1);
      tf->SetBranchAddress("AtPIDEvent", &pe);

      long before = (long)vx.size();
      Long64_t N = std::min(ts->GetEntries(), tf->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i); tf->GetEntry(i);
         double thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != pdg) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp > 0) thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thT < 0) continue;
         ++nGen;
         if (!pe || !pe->GetEntriesFast()) continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev) continue;
         double bd = 1e9, bx = 0, by = 0; bool got = false;
         for (auto &s : ev->GetSpyral()) {
            if (!s.valid) continue;
            double d = std::fabs(180.0 - s.polar * TMath::RadToDeg() - thT);
            if (d < bd) { bd = d; bx = s.sqrtdEdx; by = s.brho; got = true; }
         }
         if (got && bd < dThetaMax) { vx.push_back(bx); vy.push_back(by); }
      }
      Fs->Close(); Ff->Close();
      printf("  %-12s +%ld  (total %zu)\n", tg.Data(), (long)vx.size() - before, vx.size());
   }
   delete ta;

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   gSystem->mkdir(here + "/" + gSystem->DirName(outFile), kTRUE);
   TFile out(here + "/" + outFile, "RECREATE");
   TTree t("pts", "truth-matched simulated deuterons");
   float x, y; t.Branch("x", &x); t.Branch("y", &y);
   for (size_t i = 0; i < vx.size(); ++i) { x = vx[i]; y = vy[i]; t.Fill(); }
   t.Write(); out.Close();
   printf("\n  %zu truth-matched deuterons of %ld generated (%.1f %%)\n", vx.size(), nGen,
          100.0 * vx.size() / std::max(1L, nGen));
   printf("  wrote %s\n\n", outFile.Data());
}
