void diag49(){
  TFile* f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_reco.root");
  TTree* t=(TTree*)f->Get("cbmsim");
  TClonesArray* corr=nullptr; TClonesArray* pat=nullptr;
  t->SetBranchAddress("AtEventCorrected",&corr);
  t->SetBranchAddress("AtPatternEvent",&pat);
  t->GetEntry(49);
  AtEvent* e=(AtEvent*)corr->At(0);
  int nh=e->GetNumHits();
  // radius stats of all corrected hits
  double rmin=1e9,rmax=-1e9;
  for(int i=0;i<nh;i++){auto p=e->GetHit(i).GetPosition();double r=TMath::Hypot(p.X(),p.Y());rmin=std::min(rmin,r);rmax=std::max(rmax,r);}
  printf("AtEventCorrected: %d hits, radius %.0f..%.0f mm\n",nh,rmin,rmax);
  AtPatternEvent* pe=(AtPatternEvent*)pat->At(0);
  auto& trks=pe->GetTrackCand();
  printf("PRA tracks: %zu\n",trks.size());
  int k=0; int assigned=0;
  for(auto& tr:trks){
    int n=tr.GetHitArray().size(); assigned+=n;
    double tr_rmin=1e9,tr_rmax=-1e9;
    for(auto& hp:tr.GetHitArray()){auto p=hp->GetPosition();double r=TMath::Hypot(p.X(),p.Y());tr_rmin=std::min(tr_rmin,r);tr_rmax=std::max(tr_rmax,r);}
    printf("  track %d: %d hits, GeoTheta %.1f, GeoRadius %.1f, hit-radius %.0f..%.0f\n",
      k++,n,tr.GetGeoTheta()*TMath::RadToDeg(),tr.GetGeoRadius(),tr_rmin,tr_rmax);
  }
  printf("ASSIGNED %d / %d hits (%.0f%%), DROPPED %d (%.0f%%)\n",
    assigned,nh,100.*assigned/nh,nh-assigned,100.*(nh-assigned)/nh);
}
