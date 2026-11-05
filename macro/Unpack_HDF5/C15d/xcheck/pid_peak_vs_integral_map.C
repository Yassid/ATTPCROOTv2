// Same tracks, two charge definitions. Does a REAL pulse integral stretch the high-dE/dx end the
// way Spyral does?  No peak finding needed: define the regions on the production (peak) axis and
// see where the integral puts the same tracks.  Spyral's stretch was x1.20 at t and x1.24 at alpha.
void taumap(){
  gSystem->Load("libAtReconstruction.so");
  const char* runs[2]={"run_0016_reco.root","run_0019_reco.root"};
  AtTools::AtSpyralPID eP, eI;
  for(auto* e : {&eP,&eI}){ e->SetBField(2.85); e->SetMinPoints(30); e->SetZTieTolerance(0.0); }
  eP.SetUseTraceIntegral(false); eI.SetUseTraceIntegral(true);
  std::vector<double> vp, vi;
  for(int k=0;k<2;k++){
    TFile* F=TFile::Open(Form("/home/yassid/Data/C15d_tautest/flo2/%s",runs[k]));
    if(!F||F->IsZombie()) continue;
    TTree* t=(TTree*)F->Get("cbmsim");
    TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(!pa||pa->GetEntriesFast()==0) continue;
      auto* ev=dynamic_cast<AtPatternEvent*>(pa->At(0)); if(!ev) continue;
      for(auto& tr : ev->GetTrackCand()){
        auto rp=eP.Estimate(const_cast<AtTrack&>(tr));
        auto ri=eI.Estimate(const_cast<AtTrack&>(tr));
        if(!rp.valid||!ri.valid||rp.brho<=0||rp.sqrtdEdx<=0||ri.sqrtdEdx<=0) continue;
        if(std::hypot(rp.vertex.X(),rp.vertex.Y())>=25) continue;
        vp.push_back(rp.sqrtdEdx*rp.brho); vi.push_back(ri.sqrtdEdx*ri.brho); } }
    F->Close();
  }
  printf("paired tracks: %zu\n", vp.size());
  std::vector<double> sp=vp, si=vi;
  std::sort(sp.begin(),sp.end()); std::sort(si.begin(),si.end());
  double mp=sp[sp.size()/2], mi=si[si.size()/2];
  printf("%-14s %8s %12s %12s %9s\n","region on PEAK","n","ln PEAK","ln INTEGRAL","stretch");
  struct R{const char*n;double lo,hi;};
  R rs[4]={{"p",-0.80,-0.45},{"d",-0.12,0.12},{"t",0.28,0.52},{"alpha",0.58,0.85}};
  for(auto&r:rs){
    std::vector<double> a,b;
    for(size_t j=0;j<vp.size();++j){ double lp=log(vp[j]/mp);
      if(lp>=r.lo&&lp<r.hi){ a.push_back(lp); b.push_back(log(vi[j]/mi)); } }
    if(a.size()<50){ printf("%-14s %8zu  (too few)\n",r.n,a.size()); continue; }
    std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
    double ma=a[a.size()/2], mb=b[b.size()/2];
    printf("%-14s %8zu %12.3f %12.3f %9s\n", r.n, a.size(), ma, mb,
           fabs(ma)>0.1? Form("x%.2f",mb/ma) : "-");
  }
  printf("\n  Spyral vs ATTPCROOT for comparison:  p x0.96   d x1.00   t x1.20   alpha x1.24\n");
}
