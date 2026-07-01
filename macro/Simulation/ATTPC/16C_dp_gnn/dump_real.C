// Dump real AtEventCorrected hit clouds (for overlay canvases).
void dump_real(int nEv=600, const char* in="/mnt/f/a1975/reco_d2/run_0016_reco.root",
               const char* out="data/real_events.csv"){
  TFile*f=TFile::Open(in);TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*c=0;t->SetBranchAddress("AtEventCorrected",&c);
  FILE*fp=fopen(out,"w");fprintf(fp,"event,x,y,z,q\n");
  int n=std::min((Long64_t)nEv,t->GetEntries());
  for(int i=0;i<n;i++){t->GetEntry(i);if(!c->GetEntriesFast())continue;auto e=(AtEvent*)c->At(0);
   for(int h=0;h<e->GetNumHits();h++){auto p=e->GetHit(h).GetPosition();
    fprintf(fp,"%d,%.1f,%.1f,%.1f,%.1f\n",i,p.X(),p.Y(),p.Z(),e->GetHit(h).GetCharge());}}
  fclose(fp);printf("dumped %d real events -> %s\n",n,out);
}
