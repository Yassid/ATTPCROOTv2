/// @file gate_truth_C14.C
/// @brief Simulation analogue of pipeline/gate_events_C14.C: keep only the PRA tracks that
///        are genuinely the RECOIL PROTON, using MC truth instead of an IC + PID gate.
///
/// WHY THIS EXISTS. The data fit input is gated twice (IC on the event, Spyral PID polygon on
/// each track). The simulation had no equivalent, so every PRA track was fitted under a proton
/// hypothesis -- including the beam and the scattered 14C. The result: of 6771 fitted tracks
/// only 139 survived chi2/ndf < 5. 6605 sat at KE < 5 MeV / theta ~ 13 deg (heavy or garbage)
/// and 33 at KE > 100 MeV / theta ~ 1 deg (the beam itself). A closure test on that sample is
/// meaningless, and it is also why hyp_C14.C returned a NaN slope: it needs > 50 tracks in
/// theta_lab 55-85 deg and there were nowhere near that many real protons.
///
/// HOW. AtPulseTask::SetSaveMCInfo() (enabled in run_reco_C14.C) attaches AtHit::MCSimPoint to
/// every hit, carrying the MC A, Z and trackID. A PRA track is kept when a fraction >= purity
/// of its MC-labelled hits come from a Z=1, A=1 particle. That is a truth selection, not a
/// kinematic one, so it does NOT bias the reconstructed Ex -- which matters, because the whole
/// point of the closure test is to measure the chain's Ex bias.
///
/// Writes <outDir>/simg_reco.root with the FairRoot metadata the fitters need, exactly as
/// gate_events_C14.C does, so it can be fed to fitUKF_C14.C / fitGenfit_C14.C unchanged.
///
///   root -b -q 'gate_truth_C14.C("./data/sim_reco.root","./data/","simg",0.6,4)'

void gate_truth_C14(TString inFile = "./data/sim_reco.root", TString outDir = "./data/", TString outName = "simg",
                    Double_t purity = 0.6, Int_t minLabelled = 4)
{
   gSystem->Load("libAtReconstruction.so");

   TFile *fin = TFile::Open(inFile);
   if (!fin || fin->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", inFile.Data());
      return;
   }
   auto *t = (TTree *)fin->Get("cbmsim");
   if (!t) {
      printf("\033[1;31mno cbmsim in %s\033[0m\n", inFile.Data());
      return;
   }
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);

   TString of = outDir + outName + "_reco.root";
   TFile *fout = new TFile(of, "RECREATE", "", 1);
   TTree *nt = t->CloneTree(0);

   Long64_t nEvt = 0, nTrk = 0, nSeen = 0, nNoLabel = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!pe || pe->GetEntries() == 0)
         continue;
      auto *p = (AtPatternEvent *)pe->At(0);
      if (!p)
         continue;

      std::vector<AtTrack> protons;
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         ++nSeen;
         long nLab = 0, nProt = 0;
         for (auto &h : tr.GetHitArray()) { // HitVector holds pointers
            if (!h)
               continue;
            const auto &mc = h->GetMCSimPointArray();
            if (mc.empty())
               continue;
            // the first MC point is the dominant contributor to the hit
            ++nLab;
            if (mc[0].Z == 1 && mc[0].A == 1)
               ++nProt;
         }
         if (nLab < minLabelled) {
            ++nNoLabel;
            continue;
         }
         if (double(nProt) / double(nLab) >= purity)
            protons.push_back(tr);
      }
      if (protons.empty())
         continue;
      nTrk += protons.size();
      p->SetTrackCand(std::move(protons));
      nt->Fill();
      ++nEvt;
   }
   nt->Write();

   // FairRoot metadata, same as gate_events_C14.C, or the fitters will not read the file
   TList bl;
   for (auto b : *nt->GetListOfBranches())
      bl.Add(new TObjString(b->GetName()));
   fout->cd();
   bl.Write("BranchList", TObject::kSingleKey);
   TList etb;
   etb.Write("TimeBasedBranchList", TObject::kSingleKey);
   if (auto *fh = fin->Get("FileHeader")) {
      fout->cd();
      fh->Write("FileHeader");
   }
   if (auto *cb = fin->Get("cbmout")) {
      fout->cd();
      cb->Write("cbmout", TObject::kSingleKey);
   }
   fout->Write("", TObject::kOverwrite);

   printf("\n\033[1;33m=== truth proton gate ===\033[0m\n");
   printf("  PRA tracks seen        : %lld\n", nSeen);
   printf("  rejected, < %d MC hits : %lld\n", minLabelled, nNoLabel);
   printf("  kept as protons        : \033[1;36m%lld\033[0m  in %lld events  (purity >= %.2f)\n", nTrk, nEvt, purity);
   printf("  -> %s\n\n", of.Data());
   fout->Close();
   fin->Close();
}
