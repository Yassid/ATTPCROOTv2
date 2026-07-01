void dump_attpc_cloud(){
  gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_max.root");
  auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
  std::ofstream csv("attpc_cloud.csv"); csv<<"eid,x,y,z,q\n";
  long n=t->GetEntries();
  for(long i=0;i<n;i++){ t->GetEntry(i); auto e=(AtEvent*)ar->At(0); if(!e)continue;
    int eid=e->GetEventID();
    for(const auto&h:e->GetHits())
      csv<<eid<<","<<h->GetPosition().X()<<","<<h->GetPosition().Y()<<","<<h->GetPosition().Z()<<","<<h->GetCharge()<<"\n";
  }
  csv.close(); printf("dumped %ld ATTPCROOT events\n",n);
}
