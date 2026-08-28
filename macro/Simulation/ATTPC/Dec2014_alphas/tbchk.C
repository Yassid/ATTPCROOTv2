void tbchk(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  const int NB=8; std::vector<double> TB[NB];
  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    for(auto&h:*e->GetHitArray()){
      auto&m=h.GetMCSimPointArray(); if(m.empty()) continue;
      int pid=m[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p=(AtTpcPoint*)pts->At(pid);
      double tz=0.5*(p->GetZIn()+p->GetZOut())*10;
      int b=(int)(tz/125.); if(b<0||b>=NB) continue;
      TB[b].push_back(h.GetTimeStamp());
    }
  }
  printf("RES %8s %8s %12s %14s\n","truth z","N","median tb","expected tb");
  for(int b=0;b<NB;b++){
    if(TB[b].size()<50) continue;
    std::sort(TB[b].begin(),TB[b].end());
    double z=b*125+62;
    printf("RES %8.0f %8zu %12.0f %14.0f\n",z,TB[b].size(),TB[b][TB[b].size()/2], z/22.35/0.16);
  }
}
