// Scan all events: fraction of AtEventCorrected hits assigned to PRA tracks.
void scan_assign(){
  TFile* f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_reco.root");
  TTree* t=(TTree*)f->Get("cbmsim");
  TClonesArray* corr=nullptr; TClonesArray* pat=nullptr;
  t->SetBranchAddress("AtEventCorrected",&corr);
  t->SetBranchAddress("AtPatternEvent",&pat);
  Long64_t N=t->GetEntries();
  auto* hFrac=new TH1F("hFrac","PRA assigned fraction (events with >=50 hits);assigned/total;events",50,0,1);
  FILE* fp=fopen("/home/yassid/spyral_d2/pra_assign.csv","w");
  fprintf(fp,"event,nhits,assigned,ntracks,frac\n");
  int nsub=0,ndrop=0,nbigdrop=0;
  for(Long64_t i=0;i<N;i++){
    t->GetEntry(i);
    if(!corr||corr->GetEntriesFast()==0) continue;
    AtEvent* e=(AtEvent*)corr->At(0); int nh=e->GetNumHits();
    int assigned=0,ntr=0;
    if(pat&&pat->GetEntriesFast()>0){auto pe=(AtPatternEvent*)pat->At(0);
      for(auto& tr:pe->GetTrackCand()){assigned+=tr.GetHitArray().size();ntr++;}}
    double fr= nh>0? (double)assigned/nh : 0;
    fprintf(fp,"%lld,%d,%d,%d,%.3f\n",i,nh,assigned,ntr,fr);
    if(nh>=50){nsub++; hFrac->Fill(fr); if(fr<0.5)ndrop++; if(fr<0.3)nbigdrop++;}
    if(i%10000==0)printf("  %lld/%lld\n",i,N);
  }
  fclose(fp);
  printf("\n=== PRA assignment scan (events with >=50 corrected hits) ===\n");
  printf("substantial events: %d\n",nsub);
  printf("assigned <50%% of cloud: %d  (%.1f%%)\n",ndrop,100.*ndrop/nsub);
  printf("assigned <30%% of cloud: %d  (%.1f%%)\n",nbigdrop,100.*nbigdrop/nsub);
  printf("mean assigned fraction: %.2f, median ~ %.2f\n",hFrac->GetMean(),hFrac->GetBinCenter(hFrac->GetMaximumBin()));
  auto* c=new TCanvas("c","",700,500); hFrac->Draw("hist");
  c->SaveAs("plots/pra_assign_fraction.png");
}
