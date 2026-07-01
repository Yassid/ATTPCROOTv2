void check_noisy(){
  gSystem->Load("libAtReconstruction.so");
  std::vector<int> noisy={7777,2908,5837,566,7970,7494};
  for(TString tag : {"max","mfspyral"}){
    auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_"+tag+".root");
    auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
    long n=std::min((long)t->GetEntries(),200L);
    std::map<int,int> cnt;
    for(long i=0;i<n;i++){ t->GetEntry(i); auto e=(AtEvent*)ar->At(0); if(!e)continue;
      std::set<int> seen; for(const auto&h:e->GetHits()){int p=h->GetPadNum(); if(seen.insert(p).second)cnt[p]++;} }
    printf("%-9s (thr%s): ",tag.Data(), tag=="max"?"20":"40+prom");
    for(int p:noisy) printf("pad%d=%d%% ",p,(int)(100.0*cnt[p]/n));
    printf("\n"); f->Close();
  }
}
