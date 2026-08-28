void res2(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  const int NB=8; std::vector<double> U[NB],C[NB];
  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    for(auto&h:*e->GetHitArray()){
      auto&mcv=h.GetMCSimPointArray(); if(mcv.empty()) continue;
      int pid=mcv[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p=(AtTpcPoint*)pts->At(pid);
      double tx=0.5*(p->GetXIn()+p->GetXOut())*10, ty=0.5*(p->GetYIn()+p->GetYOut())*10;
      double tz=0.5*(p->GetZIn()+p->GetZOut())*10;    // truth z in mm = drift distance
      int b=(int)(tz/125.); if(b<0||b>=NB) continue;
      TVector3 u=h.GetPosition(), c=h.GetPositionCorr();
      U[b].push_back(TMath::Hypot(u.X()-tx,u.Y()-ty));
      C[b].push_back(TMath::Hypot(c.X()-tx,c.Y()-ty));
    }
  }
  printf("RESULT %8s %8s %10s %10s\n","z[mm]","N","uncorr","corrected");
  for(int b=0;b<NB;b++){
    if(U[b].size()<50) continue;
    std::sort(U[b].begin(),U[b].end()); std::sort(C[b].begin(),C[b].end());
    printf("RESULT %8d %8zu %10.1f %10.1f\n",b*125+62,U[b].size(),U[b][U[b].size()/2],C[b][C[b].size()/2]);
  }
}
