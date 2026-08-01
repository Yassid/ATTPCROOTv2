/// @file pid_plane_cache.C
/// @brief Cache the FULL, UNGATED Spyral PID landscape for one a1975 proton-target run.
///
/// Drawing a gate from the fit output is circular: run_XXXX_genfitter_*.root only contains tracks
/// that already passed the gate being refined, so the bands outside it are invisible. This reads
/// the PATTERN tracks in run_XXXX_reco.root instead and recomputes the PID for every one of them,
/// giving the whole proton / deuteron / triton landscape.
///
/// The recomputation mirrors AtGenfitter's in-fitter gate exactly: AtSpyralPID with SetBField(|B|)
/// and all other parameters default (AtGenfitter.cxx:121-122), so a gate drawn on this plane means
/// the same thing when the fitter later applies it.
///
/// Columns: sqrtdedx, brho, dedx, polar[deg], npoints, nclusters, run, event, trackid.
/// (p,p') and (p,d) share the same runs and the same landscape -- only the gate polygon differs.
///
///   root -b -q 'pid/pid_plane_cache.C("run_0106","/tmp/pidplane_0106.root")'

void pid_plane_cache(TString run, TString outFile, TString inDir = "/mnt/f/a1975/reco/",
                     double bField = 2.85, Long64_t maxEvents = -1)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   TString fr = inDir + run + "_reco.root";
   if (gSystem->AccessPathName(fr)) { printf("skip %s (no reco)\n", run.Data()); return; }

   TFile *Fr = TFile::Open(fr);
   TTree *tr = (TTree *)Fr->Get("cbmsim");
   tr->SetBranchStatus("*", 0);
   tr->SetBranchStatus("AtPatternEvent*", 1);
   TClonesArray *pa = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pa);

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));

   int runNo = TString(run(run.Length()-4, 4)).Atoi();
   TFile out(outFile, "RECREATE");
   float sqrtdedx, brho, dedx, polar; int npoints, nclusters, rn = runNo, evt, tid;
   TTree *pl = new TTree("pl", "ungated Spyral PID landscape");
   pl->Branch("sqrtdedx",&sqrtdedx); pl->Branch("brho",&brho); pl->Branch("dedx",&dedx);
   pl->Branch("polar",&polar); pl->Branch("npoints",&npoints); pl->Branch("nclusters",&nclusters);
   pl->Branch("run",&rn); pl->Branch("event",&evt); pl->Branch("trackid",&tid);

   Long64_t N = tr->GetEntries();
   if (maxEvents > 0 && maxEvents < N) N = maxEvents;
   long nTrk = 0, nValid = 0;
   for (Long64_t i = 0; i < N; ++i) {
      tr->GetEntry(i);
      if (pa->GetEntries() == 0) continue;
      evt = (int)i;
      for (auto &trk : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         ++nTrk;
         auto r = spy.Estimate(const_cast<AtTrack &>(trk));
         if (!r.valid) continue;
         sqrtdedx = r.sqrtdEdx; brho = r.brho; dedx = r.dEdx;
         polar = r.polar*TMath::RadToDeg(); npoints = r.nPoints;
         nclusters = r.nClusters;   // stamped by AtSpyralPID itself; no pattern lookup needed
         tid = r.trackID;
         pl->Fill(); ++nValid;
      }
   }
   out.cd(); pl->Write(); out.Close(); Fr->Close();
   printf("%s: %ld tracks -> %ld with valid PID -> %s\n", run.Data(), nTrk, nValid, outFile.Data());
}
