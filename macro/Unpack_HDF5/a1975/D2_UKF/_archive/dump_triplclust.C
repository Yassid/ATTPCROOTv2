void dump_triplclust(){ gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_multifit_reco.root");
  auto t=(TTree*)f->Get("cbmsim"); TClonesArray*pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
  std::ofstream o("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/triplclust_clusters.csv");
  o<<"event,x,y,z,cluster\n"; long n=t->GetEntries();
  for(long i=0;i<n;i++){ t->GetEntry(i); auto pe=(AtPatternEvent*)pa->At(0); if(!pe)continue;
    auto&tr=pe->GetTrackCand();
    for(size_t k=0;k<tr.size();k++)
      for(auto&h:tr[k].GetHitArray()){ auto p=h->GetPosition(); o<<i<<","<<p.X()<<","<<p.Y()<<","<<p.Z()<<","<<k<<"\n"; } }
  o.close(); printf("dumped triplclust clusters from %ld events\n",n);}
