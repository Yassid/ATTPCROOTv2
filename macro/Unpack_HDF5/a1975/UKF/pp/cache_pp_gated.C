/// @file cache_pp_gated.C
/// @brief Cache one run's (p,p') proton kinematics with the proton PID gate applied, so a UKF
///        sample can be compared like-for-like against the already-gated GENFIT production.
///
/// Why this exists: the genfit production ran the gate INSIDE the fitter
/// (AtGenfitter::SetPIDGate -> pid/proton_band.json), so run_XXXX_genfitter_pp.root already holds
/// only gated tracks. The UKF production did not -- run_XXXX_ukf.root has no AtPIDEvent branch and
/// keeps every fitted track, so caching both with pp/cache_pp_run.C compares a gated sample against
/// an ungated one (run_0106: 7772 vs 14293). Run THIS for the UKF side; genfit needs no re-gating.
///
/// The obvious shortcut -- read the per-run run_XXXX_pid.root -- does NOT work: those files predate
/// AtSpyralResult::trackID, so every entry carries trackID = -1, and the PID vector is not parallel
/// to the tracks either (event 4 of run_0106 has a fitted track of ID 2 but only one PID entry).
/// There is no way to associate them. So the PID is RECOMPUTED here from the pattern tracks in
/// run_XXXX_reco.root, where the track ID is unambiguous because we iterate the tracks ourselves.
///
/// The recomputation mirrors AtGenfitter exactly: AtSpyralPID with SetBField(|B|) and every other
/// parameter left at its default, which is all the fitter sets (AtGenfitter.cxx:121-122).
///
///   root -b -q 'pp/cache_pp_gated.C("run_0106","/tmp/g_run0106.root")'

void cache_pp_gated(TString run, TString outFile, TString inDir = "/mnt/f/a1975/reco/",
                    TString suffix = "_ukf", TString gateFile = "pid/proton_band.json",
                    double bField = 2.85, int icTbLo = 1000, int icTbHi = 1350)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   TString ft = inDir + run + suffix + ".root", fr = inDir + run + "_reco.root",
           ff = inDir + run + "_FRIB.root";
   for (auto p : {ft, fr, ff})
      if (gSystem->AccessPathName(p)) { printf("skip %s (no %s)\n", run.Data(), p.Data()); return; }

   TFile *Ft = TFile::Open(ft); TTree *tt = (TTree *)Ft->Get("cbmsim");
   TFile *Fr = TFile::Open(fr); TTree *tr = (TTree *)Fr->Get("cbmsim");
   TFile *Fc = TFile::Open(ff); TTree *tc = (TTree *)Fc->Get("cbmsim");
   // the reco file also carries the (huge) raw/event branches; only the pattern is needed
   tr->SetBranchStatus("*", 0);
   tr->SetBranchStatus("AtPatternEvent*", 1);
   TClonesArray *te = nullptr, *pa = nullptr, *re = nullptr;
   tt->SetBranchAddress("AtTrackingEvent", &te);
   tr->SetBranchAddress("AtPatternEvent", &pa);
   tc->SetBranchAddress("AtRawEvent", &re);

   auto gate = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));   // exactly what AtGenfitter sets; all else default

   int runNo = TString(run(run.Length()-4, 4)).Atoi();
   TFile out(outFile, "RECREATE");
   float ke, theta, vz, chi2ndf, ic; int rn = runNo;
   TTree *pk = new TTree("pk", "PID-gated proton kinematics");
   pk->Branch("ke",&ke); pk->Branch("theta",&theta); pk->Branch("vz",&vz);
   pk->Branch("chi2ndf",&chi2ndf); pk->Branch("ic",&ic); pk->Branch("run",&rn);

   Long64_t N = std::min({tt->GetEntries(), tr->GetEntries(), tc->GetEntries()});
   long nTot = 0, nKept = 0, nNoTrk = 0;
   for (Long64_t i = 0; i < N; ++i) {
      tc->GetEntry(i); ic = -1;
      if (re->GetEntries() > 0) {
         auto *raw = (AtRawEvent *)re->At(0);
         if (raw && !raw->GetGenTraces().empty()) {
            auto &adc = raw->GetGenTraces()[0]->GetADC();
            double mx = -1e9;
            for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b) mx = std::max(mx, adc[b]);
            ic = mx;
         }
      }
      tt->GetEntry(i);
      if (te->GetEntries() == 0) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev || ev->GetFittedTracks().empty()) continue;

      tr->GetEntry(i);
      if (pa->GetEntries() == 0) continue;
      auto &cand = ((AtPatternEvent *)pa->At(0))->GetTrackCand();

      // PID per pattern track, keyed by the SAME id the fitter stamps on AtFittedTrack
      std::map<int, bool> pass;
      for (auto &trk : cand) {
         auto r = spy.Estimate(const_cast<AtTrack &>(trk));
         pass[trk.GetTrackID()] = r.valid && gate.IsInside(r.sqrtdEdx, r.brho);
      }

      for (auto &fitted : ev->GetFittedTracks()) {
         if (!fitted) continue;
         auto &m = fitted->GetTrackMetadata();
         if (!m) continue;
         ++nTot;
         auto it = pass.find(fitted->GetTrackID());
         if (it == pass.end()) { ++nNoTrk; continue; }
         if (!it->second) continue;
         auto &k = fitted->GetKinematics();
         double ndf = m->GetNdf();
         ke = k.kineticEnergy; theta = k.theta*TMath::RadToDeg(); vz = fitted->GetVertex(0).Z();
         chi2ndf = ndf > 0 ? m->GetChi2()/ndf : 1e9;
         pk->Fill(); ++nKept;
      }
   }
   out.cd(); pk->Write(); out.Close();
   Ft->Close(); Fr->Close(); Fc->Close();
   printf("%s%s: %ld fitted -> %ld gated (%ld unmatched) -> %s\n",
          run.Data(), suffix.Data(), nTot, nKept, nNoTrk, outFile.Data());
}
