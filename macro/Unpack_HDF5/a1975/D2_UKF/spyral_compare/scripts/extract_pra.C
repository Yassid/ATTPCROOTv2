// Extract ATTPCROOT PRA clustering (AtPatternEvent) from run_0016_reco.root
// to a CSV for one-to-one comparison vs Spyral. One row per PRA track.
void extract_pra(const char* in="/mnt/f/a1975/reco_d2/run_0016_reco.root",
                 const char* out="/home/yassid/spyral_d2/pra_run0016.csv") {
   TFile* f = TFile::Open(in);
   TTree* t = (TTree*)f->Get("cbmsim");
   TClonesArray* patArr = nullptr;
   t->SetBranchAddress("AtPatternEvent", &patArr);
   Long64_t N = t->GetEntries();
   FILE* fp = fopen(out, "w");
   fprintf(fp, "event,ntracks,track_id,geo_theta,geo_phi,geo_radius,nhits,zmin,zmax,zmean\n");
   Long64_t ntrk_total = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!patArr || patArr->GetEntriesFast() == 0) continue;
      AtPatternEvent* pe = (AtPatternEvent*)patArr->At(0);
      auto& tracks = pe->GetTrackCand();
      int nt = tracks.size();
      for (int k = 0; k < nt; ++k) {
         auto& tr = tracks[k];
         auto& hits = tr.GetHitArray();
         int nh = hits.size();
         double zmin = 1e9, zmax = -1e9, zsum = 0;
         for (auto& h : hits) {
            double z = h->GetPosition().Z();
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
            zsum += z;
         }
         double zmean = nh > 0 ? zsum / nh : 0;
         fprintf(fp, "%lld,%d,%d,%.4f,%.4f,%.4f,%d,%.2f,%.2f,%.2f\n",
                 i, nt, tr.GetTrackID(),
                 tr.GetGeoTheta()*TMath::RadToDeg(),
                 tr.GetGeoPhi()*TMath::RadToDeg(),
                 tr.GetGeoRadius(), nh, zmin, zmax, zmean);
         ntrk_total++;
      }
      if (i % 5000 == 0) printf("  event %lld / %lld\n", i, N);
   }
   fclose(fp);
   printf("DONE: %lld events, %lld PRA tracks -> %s\n", N, ntrk_total, out);
}
