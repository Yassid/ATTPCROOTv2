void check_pad7777(){
  gSystem->Load("libAtReconstruction.so");
  for(TString tag : {"max","mfspyral"}){
    auto f=TFile::Open("/mnt/f/a1975/reco_d2/run_0300_psa_"+tag+".root");
    auto t=(TTree*)f->Get("cbmsim"); TClonesArray*ar=nullptr; t->SetBranchAddress("AtEventH",&ar);
    printf("%-9s pad 7777 hits per event: ", tag.Data());
    for(int i=0;i<6;i++){ t->GetEntry(i); auto e=(AtEvent*)ar->At(0);
      int nh=0; double tb=-1; for(const auto&h:e->GetHits()) if(h->GetPadNum()==7777){nh++; tb=h->GetTimeStamp();}
      printf("ev%d=%d%s ",i,nh, nh?Form("(TB%.0f)",tb):""); }
    printf("\n"); f->Close();
  }
}
