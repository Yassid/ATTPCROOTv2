void dump_mf20(){ gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_multifit.root");
  auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
  std::ofstream o("mf20_cloud.csv"); o<<"eid,x,y,z,q\n"; long n=t->GetEntries();
  for(long i=0;i<n;i++){t->GetEntry(i);auto e=(AtEvent*)ar->At(0);if(!e)continue;int id=e->GetEventID();
    for(const auto&h:e->GetHits())o<<id<<","<<h->GetPosition().X()<<","<<h->GetPosition().Y()<<","<<h->GetPosition().Z()<<","<<h->GetCharge()<<"\n";}
  o.close();printf("dumped %ld events -> mf20_cloud.csv\n",n);}
