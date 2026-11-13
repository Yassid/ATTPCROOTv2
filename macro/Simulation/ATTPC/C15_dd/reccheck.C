void reccheck(){
  gSystem->Load("libAtReconstruction.so");
  TFile f("/home/yassid/Data/C15dd_sim/test_reco.root");
  TTree* t=(TTree*)f.Get("cbmsim");
  printf("entries %lld   branches:",t->GetEntries());
  TObjArray* b=t->GetListOfBranches();
  for(int i=0;i<b->GetEntries();i++) printf(" %s",b->At(i)->GetName());
  printf("\n");
  TClonesArray* ev=nullptr; TClonesArray* pa=nullptr;
  if(t->GetBranch("AtEventH")) t->SetBranchAddress("AtEventH",&ev);
  if(t->GetBranch("AtPatternEvent")) t->SetBranchAddress("AtPatternEvent",&pa);
  long nh=0,ntr=0,nev=0,withTrk=0; double hitsPerTrk=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i); nev++;
    if(ev&&ev->GetEntriesFast()){ auto* e=(AtEvent*)ev->At(0); nh+=e->GetNumHits(); }
    if(pa&&pa->GetEntriesFast()){ auto* p=(AtPatternEvent*)pa->At(0);
      auto& tr=p->GetTrackCand(); ntr+=tr.size(); if(tr.size()) withTrk++;
      for(auto& x:tr) hitsPerTrk+=x.GetHitArray().size(); } }
  printf("  events %ld   hits/event %.1f   tracks %ld (%.0f%% of events)   hits/track %.1f\n",
         nev, nev?(double)nh/nev:0, ntr, nev?100.0*withTrk/nev:0, ntr?hitsPerTrk/ntr:0);
}
