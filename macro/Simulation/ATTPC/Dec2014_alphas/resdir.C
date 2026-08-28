// Track-DIRECTION validation of the tilt/Lorentz correction against MC truth.
//
// Hit-to-point distance has a floor set by pad granularity and track obliquity that has
// nothing to do with the drift correction. Direction removes it: fitting many hits
// averages the granularity down, and direction is the observable every physics result
// here actually rests on.
//
// Prediction fixed BEFORE looking, so the sign conventions cannot be reinterpreted after
// the fact:
//   * the UNCORRECTED direction error must grow with the track's mean drift distance,
//     because the shear is proportional to drift time and therefore tilts the fitted line;
//   * the CORRECTED direction error must be small and flat in drift.
// If the uncorrected error does not grow, the setup is wrong again and the result should
// be discarded rather than reinterpreted.
void resdir(){
  TFile*fs=TFile::Open("data/attpcsim_in.root"); TTree*ts=(TTree*)fs->Get("cbmsim");
  TClonesArray*pts=nullptr; ts->SetBranchAddress("AtTpcPoint",&pts);
  TFile*fd=TFile::Open("data/digi_driftON.root"); TTree*td=(TTree*)fd->Get("cbmsim");
  TClonesArray*ev=nullptr; td->SetBranchAddress("AtEventH",&ev);
  const double ZPAD=1000.;
  const int NB=6; std::vector<double> EU[NB],EC[NB];

  auto fitDir=[](std::vector<TVector3>&P)->TVector3{
    TVector3 c(0,0,0); for(auto&p:P) c+=p; c*= 1.0/P.size();
    double xx=0,xy=0,xz=0,yy=0,yz=0,zz=0;
    for(auto&p:P){ TVector3 d=p-c;
      xx+=d.X()*d.X(); xy+=d.X()*d.Y(); xz+=d.X()*d.Z();
      yy+=d.Y()*d.Y(); yz+=d.Y()*d.Z(); zz+=d.Z()*d.Z(); }
    TMatrixDSym M(3); M(0,0)=xx;M(0,1)=xy;M(0,2)=xz;M(1,0)=xy;M(1,1)=yy;M(1,2)=yz;M(2,0)=xz;M(2,1)=yz;M(2,2)=zz;
    TMatrixDSymEigen eig(M); TVectorD val=eig.GetEigenValues(); TMatrixD vec=eig.GetEigenVectors();
    int im=0; for(int k=1;k<3;k++) if(val[k]>val[im]) im=k;
    TVector3 d(vec(0,im),vec(1,im),vec(2,im)); return d.Unit();
  };

  for(int i=0;i<td->GetEntries();i++){
    td->GetEntry(i); ts->GetEntry(i);
    AtEvent*e=(AtEvent*)ev->At(0); if(!e) continue;
    // group hits by the MC track they came from
    std::map<int,std::vector<TVector3>> hu,hc; std::map<int,std::vector<TVector3>> tru;
    for(auto&h:*e->GetHitArray()){
      auto&m=h.GetMCSimPointArray(); if(m.empty()) continue;
      int pid=m[0].pointID; if(pid<0||pid>=pts->GetEntries()) continue;
      AtTpcPoint*p=(AtTpcPoint*)pts->At(pid);
      int trk=p->GetTrackID();
      hu[trk].push_back(h.GetPosition());
      hc[trk].push_back(h.GetPositionCorr());
      tru[trk].push_back(TVector3(0.5*(p->GetXIn()+p->GetXOut())*10,
                                  0.5*(p->GetYIn()+p->GetYOut())*10,
                                  ZPAD-0.5*(p->GetZIn()+p->GetZOut())*10));
    }
    for(auto&kv:tru){
      int trk=kv.first; if(kv.second.size()<20) continue;
      TVector3 dT=fitDir(kv.second), dU=fitDir(hu[trk]), dC=fitDir(hc[trk]);
      if(dU.Dot(dT)<0) dU=-dU; if(dC.Dot(dT)<0) dC=-dC;
      double mean=0; for(auto&p:kv.second) mean+=p.Z(); mean/=kv.second.size();
      int b=(int)(mean/170.); if(b<0||b>=NB) continue;
      EU[b].push_back(TMath::ACos(std::min(1.0,dU.Dot(dT)))*180/TMath::Pi());
      EC[b].push_back(TMath::ACos(std::min(1.0,dC.Dot(dT)))*180/TMath::Pi());
    }
  }
  printf("RES %10s %7s %14s %14s\n","drift[mm]","Ntrk","uncorr[deg]","corrected[deg]");
  for(int b=0;b<NB;b++){
    if(EU[b].size()<10) continue;
    auto med=[](std::vector<double>v){std::sort(v.begin(),v.end());return v[v.size()/2];};
    printf("RES %10.0f %7zu %14.2f %14.2f\n",b*170.+85,EU[b].size(),med(EU[b]),med(EC[b]));
  }
}
