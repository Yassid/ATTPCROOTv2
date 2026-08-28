void res3(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  const double ZPAD=1000., VXY=0.1843, VZ=2.2348;
  const int NB=8; std::vector<double> U[NB],C[NB]; std::vector<double> TBv[NB];
  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    for(auto&h:*e->GetHitArray()){
      auto&m=h.GetMCSimPointArray(); if(m.empty()) continue;
      int pid=m[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p=(AtTpcPoint*)pts->At(pid);
      double tx=0.5*(p->GetXIn()+p->GetXOut())*10, ty=0.5*(p->GetYIn()+p->GetYOut())*10;
      // DRIFT DISTANCE, i.e. distance from the pad plane -- the same convention
      // AtClusterizeTask uses (z = ZPadPlane - z_geant).
      double drift = ZPAD - 0.5*(p->GetZIn()+p->GetZOut())*10;
      int b=(int)(drift/125.); if(b<0||b>=NB) continue;
      TVector3 u=h.GetPosition(), c=h.GetPositionCorr();
      U[b].push_back(TMath::Hypot(u.X()-tx,u.Y()-ty));
      C[b].push_back(TMath::Hypot(c.X()-tx,c.Y()-ty));
      TBv[b].push_back(h.GetTimeStamp());
    }
  }
  printf("RES %9s %7s %9s %11s %11s %10s\n","drift[mm]","N","tb(med)","shear pred","uncorrected","corrected");
  for(int b=0;b<NB;b++){
    if(U[b].size()<50) continue;
    std::sort(U[b].begin(),U[b].end()); std::sort(C[b].begin(),C[b].end()); std::sort(TBv[b].begin(),TBv[b].end());
    double d=b*125+62;
    printf("RES %9.0f %7zu %9.0f %11.1f %11.1f %10.1f\n",
      d,U[b].size(),TBv[b][TBv[b].size()/2], VXY*(d/VZ)*10., U[b][U[b].size()/2], C[b][C[b].size()/2]);
  }
}
