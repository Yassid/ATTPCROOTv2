// Dump triplclust clusters (per-hit labels) from all run_030X_multifit_reco.root into one CSV
// for GNN training. Out: 16C_dp_gnn/data/attpc_all_clusters.csv (run,event,x,y,z,q,cluster).
void dump_triplclust_all()
{
   gSystem->Load("libAtReconstruction.so");
   std::ofstream o("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/attpc_all_clusters.csv");
   o << "run,event,x,y,z,q,cluster\n";
   int runs[] = {300, 301, 302, 303, 304, 305};
   for (int run : runs) {
      TString fn = Form("/mnt/f/a1975/reco_d2/run_0%d_multifit_reco.root", run);
      if (gSystem->AccessPathName(fn)) { printf("skip run_0%d (missing reco)\n", run); continue; }
      auto f = TFile::Open(fn);
      auto t = (TTree *)f->Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pa);
      long n = t->GetEntries(), pts = 0, evcl = 0;
      for (long i = 0; i < n; i++) {
         t->GetEntry(i);
         auto pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         auto &tr = pe->GetTrackCand();
         if (tr.size() > 0) evcl++;
         for (size_t k = 0; k < tr.size(); k++)
            for (auto &h : tr[k].GetHitArray()) {
               auto p = h->GetPosition();
               o << run << "," << i << "," << p.X() << "," << p.Y() << "," << p.Z() << "," << h->GetCharge() << "," << k << "\n";
               pts++;
            }
      }
      printf("run_0%d: %ld events (%ld with clusters), %ld clustered pts\n", run, n, evcl, pts);
      f->Close();
   }
   o.close();
   printf("wrote attpc_all_clusters.csv\n");
}
