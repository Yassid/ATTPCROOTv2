void find_noisy_pads(){
  gSystem->Load("libAtReconstruction.so");
  auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_max.root");
  auto t=(TTree*)f->Get("cbmsim");
  TClonesArray *ev=nullptr,*rw=nullptr;
  t->SetBranchAddress("AtEventH",&ev); t->SetBranchAddress("AtRawEvent",&rw);
  long n=std::min((long)t->GetEntries(),200L);
  // pass 1: per-pad hit frequency + position
  std::map<int,int> cnt; std::map<int,std::pair<double,double>> pos;
  for(long i=0;i<n;i++){ t->GetEntry(i); auto e=(AtEvent*)ev->At(0); if(!e)continue;
    std::set<int> seen;
    for(const auto&h:e->GetHits()){ int p=h->GetPadNum(); if(seen.insert(p).second){cnt[p]++; pos[p]={h->GetPosition().X(),h->GetPosition().Y()};} } }
  // top always-on pads
  std::vector<std::pair<int,int>> v(cnt.begin(),cnt.end());
  std::sort(v.begin(),v.end(),[](auto&a,auto&b){return a.second>b.second;});
  printf("\n=== most frequently firing pads (of %ld events) ===\n",n);
  std::vector<int> top;
  for(int k=0;k<12 && k<(int)v.size();k++){ int p=v[k].first;
    printf("  pad %5d : %3d/%ld events (%.0f%%)  at (x=%.0f, y=%.0f)\n",p,v[k].second,n,100.0*v[k].second/n,pos[p].first,pos[p].second);
    if(k<6) top.push_back(p); }
  // pass 2: dump raw traces of top-6 always-on pads for 8 events
  std::ofstream o("noisy_traces.csv"); o<<"pad,event,tb,adc\n";
  int dumped=0;
  for(long i=0;i<n && dumped<8;i++){ t->GetEntry(i); auto raw=(AtRawEvent*)rw->At(0); if(!raw)continue;
    for(const auto&pd:raw->GetPads()){ int p=pd->GetPadNum();
      if(std::find(top.begin(),top.end(),p)!=top.end()){ auto&a=pd->GetADC();
        for(int k=0;k<512;k++) o<<p<<","<<i<<","<<k<<","<<a[k]<<"\n"; } }
    dumped++; }
  o.close(); printf("\ndumped traces for top-6 pads over 8 events -> noisy_traces.csv\n");
}
