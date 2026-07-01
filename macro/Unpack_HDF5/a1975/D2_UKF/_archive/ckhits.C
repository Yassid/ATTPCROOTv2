void ckhits(){ gSystem->Load("libAtReconstruction.so");
  for(TString tag:{"max","multifit","mfspyral"}){
    auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_"+tag+".root");
    auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
    long n=std::min((long)t->GetEntries(),40L); std::vector<int> hpe;
    for(long i=0;i<n;i++){t->GetEntry(i);auto e=(AtEvent*)ar->At(0);if(e)hpe.push_back(e->GetHits().size());}
    std::sort(hpe.begin(),hpe.end());
    printf("%-10s: hits/event median %d\n",tag.Data(),hpe[hpe.size()/2]); f->Close();}}
