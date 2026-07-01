void dump_pads(){
  gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_psa_multifit.root");
  auto t=(TTree*)f->Get("cbmsim");
  TClonesArray *ev=nullptr,*rw=nullptr;
  t->SetBranchAddress("AtEventH",&ev); t->SetBranchAddress("AtRawEvent",&rw);
  t->GetEntry(57);
  auto e=(AtEvent*)ev->At(0); auto raw=(AtRawEvent*)rw->At(0);
  std::map<int,std::vector<const AtHit*>> byPad;
  for(const auto&h:e->GetHits()) byPad[h->GetPadNum()].push_back(h.get());
  // group pads by multiplicity, sort each by primary amplitude desc
  std::map<int,std::vector<int>> byMult;
  for(auto&kv:byPad) byMult[std::min((int)kv.second.size(),3)].push_back(kv.first);
  for(auto&kv:byMult) std::sort(kv.second.begin(),kv.second.end(),[&](int a,int b){
     double ma=0,mb=0; for(auto h:byPad[a])ma=std::max(ma,h->GetCharge());
     for(auto h:byPad[b])mb=std::max(mb,h->GetCharge()); return ma>mb;});
  // pick 2 each of mult 1,2,3
  std::vector<int> pads;
  for(int m=1;m<=3;m++) for(int k=0;k<2 && k<(int)byMult[m].size();k++) pads.push_back(byMult[m][k]);
  std::ofstream tr("pad_traces.csv"); tr<<"pad,mult,tb,adc\n";
  std::ofstream hi("pad_hits.csv"); hi<<"pad,peakTB,amp\n";
  for(int pn:pads){
    int mult=byPad[pn].size();
    for(const auto&p:raw->GetPads()) if(p->GetPadNum()==pn){ auto&a=p->GetADC();
      for(int k=0;k<512;k++) tr<<pn<<","<<mult<<","<<k<<","<<a[k]<<"\n"; }
    for(auto h:byPad[pn]) hi<<pn<<","<<h->GetTimeStamp()<<","<<h->GetCharge()<<"\n";
    printf("pad %d mult %d\n",pn,mult);
  }
  tr.close(); hi.close();
}
