// Is AtHit::GetTraceIntegral() an independent measurement, or a fixed multiple of GetCharge()?
// If the per-hit ratio is constant, ATTPCROOT's "integral" carries no extra information and the
// peak-vs-integral comparison was never testable inside ATTPCROOT.
void qratio(){
  gSystem->Load("libAtReconstruction.so");
  TFile* F=TFile::Open("/home/yassid/Data/C15d_reco_mf/run_0016_reco.root");
  TTree* t=(TTree*)F->Get("cbmsim");
  TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
  std::vector<double> r; std::vector<double> qs;
  for(Long64_t i=0;i<std::min(t->GetEntries(),(Long64_t)400);i++){ t->GetEntry(i);
    if(!pa||pa->GetEntriesFast()==0) continue;
    auto* ev=dynamic_cast<AtPatternEvent*>(pa->At(0)); if(!ev) continue;
    for(auto& tr : ev->GetTrackCand())
      for(auto& h : tr.GetHitArray()){ if(!h) continue;
        double q=h->GetCharge(), ti=h->GetTraceIntegral();
        if(q>0&&ti>0){ r.push_back(ti/q); qs.push_back(q); } }
  }
  std::sort(r.begin(),r.end());
  printf("hits: %zu\n", r.size());
  printf("  integral/peak ratio: min %.4f  q05 %.4f  median %.4f  q95 %.4f  max %.4f\n",
         r.front(), r[r.size()/20], r[r.size()/2], r[r.size()*19/20], r.back());
  printf("  rms/mean: ");
  double m=0; for(double v:r) m+=v; m/=r.size();
  double s=0; for(double v:r) s+=(v-m)*(v-m); s=sqrt(s/r.size());
  printf("%.4f  -> %s\n", s/m, s/m<0.02? "CONSTANT: the integral is just a rescaled peak, no new information"
                                       : "VARIES per hit: a genuine independent measurement");
  // does the ratio depend on the charge? (that is what would make the integral superlinear)
  printf("\n  ratio vs peak charge:\n");
  std::vector<std::pair<double,double>> z;
  for(size_t i=0;i<r.size();++i) z.push_back({qs[i],r[i]});
  std::sort(z.begin(),z.end());
  for(int k=0;k<5;k++){
    size_t a=z.size()*k/5, b=z.size()*(k+1)/5; double mq=0,mr=0;
    for(size_t i=a;i<b;i++){ mq+=z[i].first; mr+=z[i].second; }
    printf("    charge quintile %d: <q> %8.1f   <integral/peak> %.4f\n", k+1, mq/(b-a), mr/(b-a));
  }
}
