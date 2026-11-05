// With the peaking time floated, is the trace integral finally an INDEPENDENT measurement?
void tauratio(){
  gSystem->Load("libAtReconstruction.so");
  const char* d[2]={"/home/yassid/Data/C15d_tautest/fix/","/home/yassid/Data/C15d_tautest/flo2/"};
  const char* n[2]={"tau FIXED (production)","tau FLOATED + integral FIXED"};
  for(int k=0;k<2;k++){
    TFile* F=TFile::Open(Form("%srun_0016_reco.root",d[k]));
    if(!F||F->IsZombie()){ printf("%s: missing\n",n[k]); continue; }
    TTree* t=(TTree*)F->Get("cbmsim");
    TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
    std::vector<double> r; std::vector<std::pair<double,double>> qr;
    for(Long64_t i=0;i<std::min(t->GetEntries(),(Long64_t)400);i++){ t->GetEntry(i);
      if(!pa||pa->GetEntriesFast()==0) continue;
      auto* ev=dynamic_cast<AtPatternEvent*>(pa->At(0)); if(!ev) continue;
      for(auto& tr : ev->GetTrackCand())
        for(auto& h : tr.GetHitArray()){ if(!h) continue;
          double q=h->GetCharge(), ti=h->GetTraceIntegral();
          if(q>0&&ti>0){ r.push_back(ti/q); qr.push_back({q,ti/q}); } }
    }
    if(r.empty()){ printf("%s: no hits\n",n[k]); continue; }
    double m=0; for(double v:r)m+=v; m/=r.size();
    double s=0; for(double v:r)s+=(v-m)*(v-m); s=sqrt(s/r.size());
    std::sort(r.begin(),r.end());
    printf("\n== %-24s %zu hits\n", n[k], r.size());
    printf("   integral/peak: median %.4f  q05 %.4f  q95 %.4f  rms/mean %.4f\n",
           r[r.size()/2], r[r.size()/20], r[r.size()*19/20], s/m);
    std::sort(qr.begin(),qr.end());
    printf("   ratio vs charge:");
    for(int j=0;j<4;j++){ size_t a=qr.size()*j/4,b=qr.size()*(j+1)/4; double mq=0,mr=0;
      for(size_t i=a;i<b;i++){mq+=qr[i].first;mr+=qr[i].second;}
      printf("  <q>%.0f:%.3f", mq/(b-a), mr/(b-a)); }
    printf("\n");
    F->Close();
  }
}
