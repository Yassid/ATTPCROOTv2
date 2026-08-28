void res(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  std::vector<double> dU,dC; long n=0;
  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    for(auto&h:*e->GetHitArray()){
      auto&mcv=h.GetMCSimPointArray(); if(mcv.empty()) continue;
      int pid=mcv[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p=(AtTpcPoint*)pts->At(pid);
      double tx=0.5*(p->GetXIn()+p->GetXOut())*10, ty=0.5*(p->GetYIn()+p->GetYOut())*10;
      TVector3 u=h.GetPosition(), c=h.GetPositionCorr();
      dU.push_back(TMath::Hypot(u.X()-tx,u.Y()-ty));
      dC.push_back(TMath::Hypot(c.X()-tx,c.Y()-ty)); n++;
    }
  }
  std::sort(dU.begin(),dU.end()); std::sort(dC.begin(),dC.end());
  printf("RESULT n=%ld\n",n);
  if(n>0){
    printf("RESULT UNCORRECTED GetPosition()     median %8.2f mm\n",dU[dU.size()/2]);
    printf("RESULT CORRECTED   GetPositionCorr() median %8.2f mm\n",dC[dC.size()/2]);
    printf("RESULT improvement %.1fx\n",dU[dU.size()/2]/dC[dC.size()/2]);
  }
}
