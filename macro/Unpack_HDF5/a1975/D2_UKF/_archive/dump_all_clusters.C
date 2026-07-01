// Dump (run,event,x,y,z,q,cluster) from a mover-join reco: full cloud (AtEventH) + per-hit track
// membership (AtPatternEvent). cluster = -1 for noise (LOF-removed + HDBSCAN noise). For GNN labels.
//   root -l -b -q 'dump_all_clusters.C("reco.root",300,"out.csv",false)'   // false = write header
//   root -l -b -q 'dump_all_clusters.C("reco.root",301,"out.csv",true)'    // true  = append
void dump_all_clusters(TString recoFile, int run, TString outcsv, bool append)
{
   gSystem->Load("libAtReconstruction.so");
   auto f = TFile::Open(recoFile);
   if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", recoFile.Data()); return; }
   auto t = (TTree *)f->Get("cbmsim");
   TClonesArray *ea = nullptr, *pa = nullptr;
   t->SetBranchAddress("AtEventH", &ea);
   t->SetBranchAddress("AtPatternEvent", &pa);
   std::ofstream o;
   if (append)
      o.open(outcsv, std::ios::app);
   else {
      o.open(outcsv);
      o << "run,event,x,y,z,q,cluster\n";
   }
   long nev = t->GetEntries();
   for (long i = 0; i < nev; i++) {
      t->GetEntry(i);
      auto ev = (AtEvent *)ea->At(0);
      if (!ev) continue;
      int nh = ev->GetNumHits();
      if (nh < 10) continue;
      auto pe = (AtPatternEvent *)pa->At(0);
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
         auto &h = ev->GetHit(j);
         auto p = h.GetPosition();
         o << run << "," << i << "," << p.X() << "," << p.Y() << "," << p.Z() << "," << h.GetCharge() << ","
           << lab[j] << "\n";
      }
   }
   o.close();
   printf("dumped run %d (%ld ev) -> %s\n", run, nev, outcsv.Data());
}
