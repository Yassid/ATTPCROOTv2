// Dump full PSA cloud + triplclust labels for a run -> 16C_dp_gnn/data/cloud_<run>.csv
// (event,x,y,z,q,cluster ; cluster=-1 means not clustered by triplclust). For the track viewer.
void dump_cloud_clusters(TString run = "run_0305", int maxEv = 400)
{
   gSystem->Load("libAtReconstruction.so");
   TString fn = "/mnt/f/a1975/reco_d2/" + run + "_multifit_reco.root";
   auto f = TFile::Open(fn);
   auto t = (TTree *)f->Get("cbmsim");
   TClonesArray *ea = nullptr, *pa = nullptr;
   t->SetBranchAddress("AtEventH", &ea);
   t->SetBranchAddress("AtPatternEvent", &pa);
   TString out = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/cloud_" + run + ".csv";
   std::ofstream o(out.Data());
   o << "event,x,y,z,q,cluster\n";
   long n = std::min((Long64_t)maxEv, t->GetEntries());
   for (long i = 0; i < n; i++) {
      t->GetEntry(i);
      auto ev = (AtEvent *)ea->At(0);
      if (!ev) continue;
      int nh = ev->GetNumHits();
      std::vector<int> lab(nh, -1);
      auto pe = (AtPatternEvent *)pa->At(0);
      if (pe) {
         auto &tr = pe->GetTrackCand();
         for (size_t k = 0; k < tr.size(); k++)
            for (auto &h : tr[k].GetHitArray()) {
               int id = h->GetHitID();
               if (id >= 0 && id < nh) lab[id] = (int)k;
            }
      }
      for (int j = 0; j < nh; j++) {
         auto &h = ev->GetHit(j);
         auto p = h.GetPosition();
         o << i << "," << p.X() << "," << p.Y() << "," << p.Z() << "," << h.GetCharge() << "," << lab[j] << "\n";
      }
   }
   o.close();
   printf("wrote %s (%ld events)\n", out.Data(), n);
}
