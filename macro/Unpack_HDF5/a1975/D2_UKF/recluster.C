// Re-cluster the AtEventH cloud of a reco file with triplclust or HDBSCAN; dump per-hit labels.
//   root -l -b -q 'recluster.C("in_reco.root","out.csv","tc",0.03,4.0)'      // triplclust a,t
//   root -l -b -q 'recluster.C("in_reco.root","out.csv","hdbscan",0,0,15,6)' // HDBSCAN mcs,ms
void recluster(TString infile, TString outcsv, TString algo, double a = 0.03, double t = 4.0, int mcs = 15,
               int ms = 6, int maxEv = 100000)
{
   gSystem->Load("libAtReconstruction.so");
   auto f = TFile::Open(infile);
   auto tr = (TTree *)f->Get("cbmsim");
   TClonesArray *ea = nullptr;
   tr->SetBranchAddress("AtEventH", &ea);
   std::ofstream o(outcsv);
   o << "event,x,y,z,label\n";
   long n = std::min((Long64_t)maxEv, tr->GetEntries());
   for (long i = 0; i < n; i++) {
      tr->GetEntry(i);
      auto ev = (AtEvent *)ea->At(0);
      if (!ev) continue;
      int nh = ev->GetNumHits();
      if (nh < 10) continue;
      std::unique_ptr<AtPatternEvent> pe;
      if (algo == "hdbscan") {
         AtPATTERN::AtTrackFinderHDBSCAN hd;
         hd.SetMinClusterSize(mcs);
         hd.SetMinSamples(ms);
         pe = hd.FindTracks(*ev);
      } else {
         AtPATTERN::AtTrackFinderTC tc;
         tc.SetClusterRadius(15.0); tc.SetClusterDistance(7.5); // REQUIRED: else tracks dropped (NaN GeoRadius)
         tc.SetScluster(0.3); tc.SetKtriplet(19); tc.SetNtriplet(2); tc.SetMcluster(15);
         tc.SetRsmooth(2); tc.SetAtriplet(a); tc.SetTcluster(t);
         try { pe = tc.FindTracks(*ev); } catch (...) {}
      }
      std::vector<int> lab(nh, -1);
      if (pe) {
         auto &trk = pe->GetTrackCand();
         for (size_t k = 0; k < trk.size(); k++)
            for (auto &h : trk[k].GetHitArray()) {
               int id = h->GetHitID();
               if (id >= 0 && id < nh) lab[id] = (int)k;
            }
      }
      for (int j = 0; j < nh; j++) {
         auto p = ev->GetHit(j).GetPosition();
         o << i << "," << p.X() << "," << p.Y() << "," << p.Z() << "," << lab[j] << "\n";
      }
   }
   o.close();
   printf("wrote %s\n", outcsv.Data());
}
