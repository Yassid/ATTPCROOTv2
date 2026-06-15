/// Dump the Spyral PID plane (sqrtdEdx, brho) of a proton-hyp genfit output, matched
/// to per-track fit quality, for auto-deriving the D2 proton gate.
///   root -b -q 'pid/dumpPID_d2.C("/tmp/run_0016_genfitter_p_pid.root","/tmp/pid_d2.txt")'
/// then: python pid/derive_proton_gate_d2.py /tmp/pid_d2.txt
void dumpPID_d2(TString in = "/tmp/run_0016_genfitter_p_pid.root", TString out = "/tmp/pid_d2.txt")
{
   gSystem->Load("libAtReconstruction.so");
   TFile *fi = TFile::Open(in);
   TTree *t = (TTree *)fi->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   t->SetBranchAddress("AtPIDEvent", &pe);
   FILE *o = fopen(out.Data(), "w");
   fprintf(o, "sqrtdEdx brho chi2ndf KE vertexR polar_deg direction\n");
   long nval = 0;
   for (Long64_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (pe->GetEntries() == 0)
         continue;
      auto *pidev = (AtPIDEvent *)pe->At(0);
      if (!pidev)
         continue;
      std::map<int, AtFittedTrack *> fmap; // fitted tracks by trackID
      if (te->GetEntries() > 0) {
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (ev)
            for (auto &ft : ev->GetFittedTracks())
               if (ft)
                  fmap[ft->GetTrackID()] = ft.get();
      }
      for (auto &sr : pidev->GetSpyral()) {
         if (!sr.valid)
            continue;
         nval++;
         double chi2ndf = -1, ke = -1, vr = -1;
         auto it = fmap.find(sr.trackID);
         if (it != fmap.end()) {
            auto *ft = it->second;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            chi2ndf = ndf > 0 ? chi2 / ndf : -1;
            ke = k.kineticEnergy;
            auto v = ft->GetVertex();
            vr = TMath::Sqrt(v.X() * v.X() + v.Y() * v.Y());
         }
         fprintf(o, "%.4f %.5f %.4f %.4f %.2f %.2f %d\n", sr.sqrtdEdx, sr.brho, chi2ndf, ke, vr,
                 sr.polar * TMath::RadToDeg(), sr.direction);
      }
   }
   fclose(o);
   printf("RESULT: dumped %ld valid PID results -> %s\n", nval, out.Data());
}
