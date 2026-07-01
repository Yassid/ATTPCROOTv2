void dump_padxy(){
  gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_max.root");
  auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
  std::map<int,std::pair<double,double>> pad;
  for(long i=0;i<std::min((long)t->GetEntries(),200L);i++){ t->GetEntry(i); auto e=(AtEvent*)ar->At(0); if(!e)continue;
    for(const auto&h:e->GetHits()) pad[h->GetPadNum()]={h->GetPosition().X(),h->GetPosition().Y()}; }
  std::ofstream o("attpc_padxy.csv"); o<<"pad,x,y\n";
  for(auto&kv:pad) o<<kv.first<<","<<kv.second.first<<","<<kv.second.second<<"\n";
  printf("dumped %zu ATTPCROOT pads\n",pad.size());
}
