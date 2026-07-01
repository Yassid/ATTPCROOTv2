void ck(){ gSystem->Load("libAtReconstruction.so");
  std::vector<int> noisy={7777,2908,5837,566};
  for(TString tag:{"max","multifit","mfspyral"}){
    auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_"+tag+".root");
    auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
    long n=std::min((long)t->GetEntries(),40L); std::map<int,int> c;
    for(long i=0;i<n;i++){t->GetEntry(i);auto e=(AtEvent*)ar->At(0);if(!e)continue;std::set<int>s;
      for(const auto&h:e->GetHits())if(s.insert(h->GetPadNum()).second)c[h->GetPadNum()]++;}
    printf("%-10s: ",tag.Data()); for(int p:noisy)printf("pad%d=%d%% ",p,(int)(100.0*c[p]/n)); printf("\n");f->Close();}}
