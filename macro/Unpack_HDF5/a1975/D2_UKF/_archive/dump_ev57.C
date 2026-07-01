void dump_ev57(){
  gSystem->Load("libAtReconstruction.so");
  std::ofstream csv("ev57_clouds.csv"); csv<<"psa,x,y,z,pad,mult,rank\n";
  for(TString tag : {"max","multifit"}){
    auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_psa_"+tag+".root");
    auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
    t->GetEntry(57); auto ev=(AtEvent*)ar->At(0);
    std::map<int,std::vector<const AtHit*>> byPad;
    for(const auto&h:ev->GetHits()) byPad[h->GetPadNum()].push_back(h.get());
    long nh=0,nmp=0;
    for(auto&kv:byPad){ auto hs=kv.second;
      std::sort(hs.begin(),hs.end(),[](auto a,auto b){return a->GetCharge()>b->GetCharge();});
      if(hs.size()>1)nmp++;
      for(size_t r=0;r<hs.size();++r){auto h=hs[r];nh++;
        csv<<tag<<","<<h->GetPosition().X()<<","<<h->GetPosition().Y()<<","<<h->GetPosition().Z()<<","
           <<h->GetPadNum()<<","<<hs.size()<<","<<r<<"\n";}}
    printf("%-9s event57: %ld hits, %lu pads, %ld multi-hit pads\n",tag.Data(),nh,byPad.size(),nmp);
    f->Close();
  }
  csv.close();
}
