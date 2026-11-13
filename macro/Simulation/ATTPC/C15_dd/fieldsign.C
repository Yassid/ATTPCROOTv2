// FIELD SIGN CHECK.  The simulation generates at -28.5 kG; the C15d fit is given +2.85 T.  If the
// sense of rotation disagrees, anything inferring direction from curvature rejects most simulated
// tracks and the acceptance measures the sign error instead of the detector.
//
// Observable: order each track's hits by z and sum the signed area  (x_i*y_{i+1} - x_{i+1}*y_i).
// Its sign is the rotation sense seen from +z, independent of where the track starts.
void sense(const char* f,const char* tag,Long64_t maxEv){
  TFile fh(f); TTree* t=(TTree*)fh.Get("cbmsim"); if(!t){printf("  %-30s no tree\n",tag);return;}
  TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
  Long64_t n=std::min(maxEv,t->GetEntries());
  long pos=0,neg=0; double sA=0;
  for(Long64_t i=0;i<n;i++){ t->GetEntry(i);
    if(!pa->GetEntriesFast()) continue;
    for(auto& tr : ((AtPatternEvent*)pa->At(0))->GetTrackCand()){
      auto hv=tr.GetHitArrayObject(); if(hv.size()<8) continue;
      std::sort(hv.begin(),hv.end(),[](const AtHit&a,const AtHit&b){
        return a.GetPosition().Z()<b.GetPosition().Z(); });
      double A=0;
      for(size_t k=0;k+1<hv.size();k++){
        auto p1=hv[k].GetPosition(), p2=hv[k+1].GetPosition();
        A += p1.X()*p2.Y() - p2.X()*p1.Y(); }
      if(A>0) pos++; else neg++; sA+=A; } }
  long tot=pos+neg;
  printf("  %-30s tracks %5ld   clockwise %5.1f %%   anticlockwise %5.1f %%   <signed area> %+.3g\n",
         tag,tot,tot?100.0*neg/tot:0,tot?100.0*pos/tot:0,tot?sA/tot:0);
}
void fieldsign(){
  gSystem->Load("libAtReconstruction.so");
  sense("/home/yassid/Data/C15d_reco_mfi/run_0017_reco.root","DATA run_0017",4000);
  sense("/home/yassid/Data/C15d_reco_mfi/run_0021_reco.root","DATA run_0021",4000);
  sense("/home/yassid/Data/C15dd_sim/g150000_reco.root","SIM  (field -28.5 kG)",150);
}
