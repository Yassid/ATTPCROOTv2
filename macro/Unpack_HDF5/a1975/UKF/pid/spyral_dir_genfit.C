/// @file spyral_dir_genfit.C
/// @brief Attach the production genfit (p,p) fit to every track of spyral_dir_audit.C's tree.
///
/// Question: where AtPRA's GeoTheta and AtSpyralPID's direction disagree, and Spyral's call puts the
/// track on the 16C(p,p) elastic locus (941 tracks, runs 0106-0108), what did genfit -- seeded from
/// AtPRA -- make of it? Match is (entry, trackID) in <run>_genfitter_pphand.root, entries aligned with
/// <run>_reco.root (30094 = 30094 on run_0106).
/// DUPLICATE-ENTRY GUARD (reference_genfitter_tree_duplicate_entries): a copied entry has an EMPTY
/// AtPIDEvent, so an entry is used only if its PID event carries a Spyral result with that trackID.
/// Kinematics: GetKinematics() (RAW, genfit slot), vertex Z as stored, as cache_pp_run.C does.
///
///   root -l -b -q 'pid/spyral_dir_genfit.C'
void spyral_dir_genfit(TString audit = "pid/plots/spyral_dir_audit.root",
                       TString out = "pid/plots/spyral_dir_genfit.root", TString inDir = "/mnt/f/a1975/reco/",
                       TString suffix = "_genfitter_pphand")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kFatal;
   TFile fa(audit);
   auto *d = (TTree *)fa.Get("d");
   int run, tid;
   Long64_t ev;
   d->SetBranchAddress("run", &run); d->SetBranchAddress("ev", &ev); d->SetBranchAddress("tid", &tid);

   TFile fo(out, "RECREATE");
   TTree *g = d->CloneTree(0);
   g->SetName("g");
   int gfit, gpid;
   double gke, gth, gvz, gchi;
   g->Branch("gfit", &gfit); g->Branch("gpid", &gpid); g->Branch("gke", &gke); g->Branch("gth", &gth);
   g->Branch("gvz", &gvz); g->Branch("gchi", &gchi);

   int curRun = -1;
   TFile *fg = nullptr;
   TTree *tg = nullptr;
   TClonesArray *te = nullptr, *tp = nullptr;
   Long64_t curEv = -1;
   for (Long64_t i = 0; i < d->GetEntries(); ++i) {
      d->GetEntry(i);
      if (run != curRun) {
         if (fg) fg->Close();
         fg = TFile::Open(Form("%s/run_%04d%s.root", inDir.Data(), run, suffix.Data()));
         tg = (TTree *)fg->Get("cbmsim");
         te = tp = nullptr;
         tg->SetBranchAddress("AtTrackingEvent", &te);
         tg->SetBranchAddress("AtPIDEvent", &tp);
         curRun = run; curEv = -1;
      }
      if (ev != curEv) { tg->GetEntry(ev); curEv = ev; }
      gfit = gpid = 0; gke = gth = gvz = gchi = -999;
      auto *pe = (tp && tp->GetEntries()) ? (AtPIDEvent *)tp->At(0) : nullptr;
      if (pe)
         for (auto &s : pe->GetSpyral()) if (s.trackID == tid) gpid = 1;
      auto *tev = (te && te->GetEntries()) ? (AtTrackingEvent *)te->At(0) : nullptr;
      if (gpid && tev)
         for (auto &ft : tev->GetFittedTracks()) {
            if (!ft || ft->GetTrackID() != tid) continue;
            auto &m = ft->GetTrackMetadata();
            if (!m) continue;
            auto &k = ft->GetKinematics();
            gfit = 1;
            gke = k.kineticEnergy; gth = k.theta * TMath::RadToDeg(); gvz = ft->GetVertex(0).Z();
            gchi = m->GetNdf() > 0 ? m->GetChi2() / m->GetNdf() : 1e9;
         }
      g->Fill();
   }
   fo.cd(); g->Write(); fo.Close();
}
