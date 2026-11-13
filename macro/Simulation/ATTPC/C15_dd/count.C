void count(){
  void* d=gSystem->OpenDirectory("/home/yassid/Data/C15dd_sim");
  std::vector<TString> f; const char* e;
  while((e=gSystem->GetDirEntry(d))){ TString s(e); if(s.BeginsWith("ex")&&s.EndsWith("_sim.root")) f.push_back(s); }
  gSystem->FreeDirectory(d); std::sort(f.begin(),f.end());
  for(auto&s:f){ TFile fh("/home/yassid/Data/C15dd_sim/"+s);
    TTree* t=(TTree*)fh.Get("cbmsim");
    printf("  %-22s %8lld events %s\n",s.Data(),t?t->GetEntries():-1,
           (t&&t->GetEntries()==20000)?"":"  <-- INCOMPLETE"); }
}
