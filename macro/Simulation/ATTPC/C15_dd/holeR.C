// Measure the a2091 beam-hole radius the way AtBeamHole.h prescribes: histogram the HIT RADIUS
// in the data and in the simulation, normalise each to its own total, and take the radius beyond
// which the sim/data ratio goes flat.  AtBeamHole.h records that 15C has never been checked.
// The C15d data reco persists only AtPatternEvent and AtPIDEvent, so compare the hits ON TRACKS
// in both -- which is the fairer comparison anyway, since those are the hits the analysis uses.
void fill(const char* f,TH1D* h,Long64_t maxEv){
  TFile fh(f); TTree* t=(TTree*)fh.Get("cbmsim"); if(!t){printf("  %s: no tree\n",f);return;}
  TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
  Long64_t n=std::min(maxEv,t->GetEntries());
  for(Long64_t i=0;i<n;i++){ t->GetEntry(i);
    if(!pa->GetEntriesFast()) continue;
    auto* pe=(AtPatternEvent*)pa->At(0);
    for(auto& tr : pe->GetTrackCand())
      for(auto& hit : tr.GetHitArrayObject()){   // HitVector holds pointers; this gives objects
        auto p=hit.GetPosition(); h->Fill(std::sqrt(p.X()*p.X()+p.Y()*p.Y())); } }
}
void holeR(){
  gSystem->Load("libAtReconstruction.so");
  TH1D hd("hd","",60,0,120), hs("hs","",60,0,120);
  fill("/home/yassid/Data/C15d_reco_mfi/run_0017_reco.root",&hd,8000);
  fill("/home/yassid/Data/C15dd_sim/hr0_reco.root",&hs,100);
  double sd=hd.Integral(),ss=hs.Integral();
  printf("  data hits %.0f   sim hits %.0f\n",sd,ss);
  if(sd<=0||ss<=0){printf("  EMPTY\n");return;}
  hd.Scale(1/sd); hs.Scale(1/ss);
  printf("   r (mm)    data frac   sim frac    sim/data\n");
  for(int b=1;b<=30;b++){
    double d=hd.GetBinContent(b), s=hs.GetBinContent(b);
    printf("   %3.0f-%3.0f   %9.5f  %9.5f  %9.2f%s\n",
           hd.GetBinLowEdge(b),hd.GetBinLowEdge(b)+hd.GetBinWidth(b),d,s,d>0?s/d:-1,
           (d>0&&s/d>3)?"   <-- sim excess":""); }
}
