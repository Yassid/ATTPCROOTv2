// Test the matching hypothesis: is the middle-drift residual an artefact of labelling a
// hit by MCSimPoint[0], when several MC points contribute to the same pad?
void res4(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  const double ZPAD=1000.;
  const int NB=8;
  std::vector<double> U0[NB],C0[NB],Uc[NB],Cc[NB]; std::vector<double> NP[NB];
  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    for(auto&h:*e->GetHitArray()){
      auto&m=h.GetMCSimPointArray(); if(m.empty()) continue;
      // (a) first point only, as before
      int pid=m[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p0=(AtTpcPoint*)pts->At(pid);
      double x0=0.5*(p0->GetXIn()+p0->GetXOut())*10, y0=0.5*(p0->GetYIn()+p0->GetYOut())*10;
      double d0=ZPAD-0.5*(p0->GetZIn()+p0->GetZOut())*10;
      // (b) energy-loss-weighted centroid over ALL contributing points
      double sx=0,sy=0,sd=0,sw=0;
      for(auto&mp:m){
        if(mp.pointID<0||mp.pointID>=pts->GetEntries()) continue;
        AtTpcPoint*p=(AtTpcPoint*)pts->At(mp.pointID);
        double w=mp.eloss>0?mp.eloss:1.0;
        sx+=w*0.5*(p->GetXIn()+p->GetXOut())*10; sy+=w*0.5*(p->GetYIn()+p->GetYOut())*10;
        sd+=w*(ZPAD-0.5*(p->GetZIn()+p->GetZOut())*10); sw+=w;
      }
      if(sw<=0) continue;
      double xc=sx/sw, yc=sy/sw, dc=sd/sw;
      int b=(int)(dc/125.); if(b<0||b>=NB) continue;
      TVector3 u=h.GetPosition(), c=h.GetPositionCorr();
      U0[b].push_back(TMath::Hypot(u.X()-x0,u.Y()-y0));
      C0[b].push_back(TMath::Hypot(c.X()-x0,c.Y()-y0));
      Uc[b].push_back(TMath::Hypot(u.X()-xc,u.Y()-yc));
      Cc[b].push_back(TMath::Hypot(c.X()-xc,c.Y()-yc));
      NP[b].push_back(m.size());
    }
  }
  const double VXY=0.1843,VZ=2.2348;
  printf("RES %9s %7s %8s %9s | %9s %9s | %9s %9s\n",
     "drift[mm]","N","<pts/hit>","shear","uncorr[0]","corr[0]","uncorr[c]","corr[c]");
  for(int b=0;b<NB;b++){
    if(U0[b].size()<50) continue;
    auto med=[](std::vector<double>v){std::sort(v.begin(),v.end());return v[v.size()/2];};
    double d=b*125+62;
    double np=0; for(double v:NP[b]) np+=v; np/=NP[b].size();
    printf("RES %9.0f %7zu %8.1f %9.1f | %9.1f %9.1f | %9.1f %9.1f\n",
      d,U0[b].size(),np,VXY*d/VZ,med(U0[b]),med(C0[b]),med(Uc[b]),med(Cc[b]));
  }
}
