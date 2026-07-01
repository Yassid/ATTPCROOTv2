// dump real AtEventH full PSA cloud for distribution comparison
void dumpRealH(int nEv=3000, const char* in="/mnt/f/a1975/reco_d2/run_0016_multifit_reco.root",
               const char* out="/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad/real_events.csv"){
  TFile*f=TFile::Open(in);TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*c=0;t->SetBranchAddress("AtEventH",&c);
  FILE*fp=fopen(out,"w");fprintf(fp,"event,x,y,z,q\n");
  int n=std::min((Long64_t)nEv,t->GetEntries());long nh=0;
  for(int i=0;i<n;i++){t->GetEntry(i);if(!c->GetEntriesFast())continue;auto e=(AtEvent*)c->At(0);
   for(int h=0;h<e->GetNumHits();h++){auto p=e->GetHit(h).GetPosition();
    fprintf(fp,"%d,%.1f,%.1f,%.1f,%.1f\n",i,p.X(),p.Y(),p.Z(),e->GetHit(h).GetCharge());nh++;}}
  fclose(fp);printf("dumped %d real events, %ld hits -> %s\n",n,nh,out);
}
