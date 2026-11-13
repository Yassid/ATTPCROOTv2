// The gain validation run_reco_C15dd.C's header demands: compare hits-per-track and
// tracks-per-event between the simulation and the DATA production it is meant to correct.
// Too low a gain silently deletes short, low-dE/dx tracks -- precisely what an acceptance measures.
void stat(const char* f,const char* tag,Long64_t maxEv){
  TFile fh(f); TTree* t=(TTree*)fh.Get("cbmsim"); if(!t){printf("  %-28s no tree\n",tag);return;}
  TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
  Long64_t n=std::min(maxEv,t->GetEntries());
  long ntr=0,wt=0; double hpt=0; std::vector<double> v;
  for(Long64_t i=0;i<n;i++){ t->GetEntry(i);
    if(!pa->GetEntriesFast()) continue;
    auto* pe=(AtPatternEvent*)pa->At(0);
    auto& tc=pe->GetTrackCand(); ntr+=tc.size(); if(tc.size()) wt++;
    for(auto& x:tc){ double h=x.GetHitArray().size(); hpt+=h; v.push_back(h); } }
  std::sort(v.begin(),v.end());
  printf("  %-28s events %6lld  with a track %5.1f %%  tracks/event %5.2f  hits/track mean %6.1f median %6.1f\n",
         tag,n,100.0*wt/n,(double)ntr/n,ntr?hpt/ntr:0,v.empty()?0:v[v.size()/2]);
}
void gainval(){
  gSystem->Load("libAtReconstruction.so");
  stat("/home/yassid/Data/C15d_reco_mfi/run_0017_reco.root","DATA run_0017 (mfi)",8000);
  stat("/home/yassid/Data/C15dd_sim/hr20_reco.root","SIM gain 1.0e4",150);
  for(const char* g : {"50000","150000","400000","1000000"})
    stat(Form("/home/yassid/Data/C15dd_sim/g%s_reco.root",g),Form("SIM gain %s",g),150);
}
