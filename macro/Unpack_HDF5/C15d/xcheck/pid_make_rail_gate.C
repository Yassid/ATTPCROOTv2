#include "ridge.C"
// Build a PROPOSED deuteron gate that follows the measured deuteron rail +-2.5 sigma, using the
// measured relative width. Written to a NEW file -- the user's own gate is not touched.
void mkgate(double nsig=2.5){
  const double K0=0.4242, rel=0.068;
  std::vector<std::pair<double,double>> up, dn;
  for(double br=0.20; br<=1.30; br+=0.02){
    double mu=K0*sqrt(betheRel(br,2,1)); if(mu<0.5||mu>44) continue;
    up.push_back({mu*(1+nsig*rel), br});
    dn.push_back({mu*(1-nsig*rel), br});
  }
  FILE* f=fopen("/home/yassid/Data/C15d_pidstudy/deuteron_rail_proposal.json","w");
  fprintf(f,"{\n    \"name\": \"deuteron_rail\",\n    \"Z\": 1,\n    \"A\": 2,\n    \"vertices\": [\n");
  // up the right edge, back down the left edge
  bool first=true;
  for(auto it=up.begin(); it!=up.end(); ++it){
    fprintf(f,"%s        [%.4f, %.4f]", first?"":",\n", it->first, it->second); first=false; }
  for(auto it=dn.rbegin(); it!=dn.rend(); ++it)
    fprintf(f,",\n        [%.4f, %.4f]", it->first, it->second);
  fprintf(f,",\n        [%.4f, %.4f]\n    ]\n}\n", up.front().first, up.front().second);
  fclose(f);
  printf("wrote deuteron_rail_proposal.json  (%zu+%zu vertices, +-%.1f sigma)\n",up.size(),dn.size(),nsig);
  // evaluate it
  TCutG cut("g",up.size()+dn.size()+1); int k=0;
  for(auto&q:up) cut.SetPoint(k++,q.first,q.second);
  for(auto it=dn.rbegin(); it!=dn.rend(); ++it) cut.SetPoint(k++,it->first,it->second);
  cut.SetPoint(k,up.front().first,up.front().second);
  TFile ff("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root");
  TTree* t=(TTree*)ff.Get("pts");
  float x,y,ic;int np;
  t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
  t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
  const double W=0.068*2.0;
  long n=0,d=0,tt=0,a=0,none=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    if(ic<1045||ic>1225||np!=1||x<=0||y<=0) continue;
    if(!cut.IsInside(x,y)) continue; n++;
    bool dD=fabs(x-K0*sqrt(betheRel(y,2,1)))/(K0*sqrt(betheRel(y,2,1)))<W;
    bool dT=fabs(x-K0*sqrt(betheRel(y,3,1)))/(K0*sqrt(betheRel(y,3,1)))<W;
    bool dA=fabs(x-K0*sqrt(betheRel(y,4,2)))/(K0*sqrt(betheRel(y,4,2)))<W;
    if(dD) d++; else if(dT) tt++; else if(dA) a++; else none++; }
  printf("  proposal selects %ld tracks:  deuteron %.1f %%  triton %.1f %%  alpha %.1f %%  bkg %.1f %%\n",
         n,100.*d/n,100.*tt/n,100.*a/n,100.*none/n);
  printf("  (saved gate selects 450931 tracks: deuteron 37.7 %%, triton 42.8 %%)\n");
}
