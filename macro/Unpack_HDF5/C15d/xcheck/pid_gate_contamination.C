#include "ridge.C"
// Refined contamination estimate: "nearest rail" assigns every track to something. Instead count
// tracks within +-2 sigma of each rail using the MEASURED relative width (sigma/mu = 0.068), and
// report the ambiguous and unassigned populations honestly.
void leak2(){
  std::vector<double> gx,gy;
  std::ifstream in("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/deuteron_C15d.json");
  std::string all((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
  size_t p=all.find("\"vertices\"");
  while(true){ size_t a=all.find('[',p+1); if(a==std::string::npos)break;
    size_t b=all.find(']',a); if(b==std::string::npos)break;
    double u,v; if(sscanf(all.substr(a+1,b-a-1).c_str()," %lf , %lf",&u,&v)==2){gx.push_back(u);gy.push_back(v);} p=b; }
  TCutG cut("d",gx.size()); for(size_t i=0;i<gx.size();++i) cut.SetPoint(i,gx[i],gy[i]);
  TFile f("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root");
  TTree* t=(TTree*)f.Get("pts");
  float x,y,ic;int np;
  t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
  t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
  const double K0=0.4242, W=0.068*2.0;   // +-2 sigma using the measured relative width
  long nTot=0,onlyD=0,onlyT=0,onlyA=0,onlyP=0,amb=0,none=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    if(ic<1045||ic>1225||np!=1||x<=0||y<=0) continue;
    if(!cut.IsInside(x,y)) continue;
    nTot++;
    bool dP=fabs(x-K0*sqrt(betheRel(y,1,1)))/(K0*sqrt(betheRel(y,1,1)))<W;
    bool dD=fabs(x-K0*sqrt(betheRel(y,2,1)))/(K0*sqrt(betheRel(y,2,1)))<W;
    bool dT=fabs(x-K0*sqrt(betheRel(y,3,1)))/(K0*sqrt(betheRel(y,3,1)))<W;
    bool dA=fabs(x-K0*sqrt(betheRel(y,4,2)))/(K0*sqrt(betheRel(y,4,2)))<W;
    int n=dP+dD+dT+dA;
    if(n==0) none++;
    else if(n>1) amb++;
    else if(dD) onlyD++; else if(dT) onlyT++; else if(dA) onlyA++; else onlyP++;
  }
  printf("inside the saved deuteron gate (IC 1045-1225, single pulse): %ld tracks\n",nTot);
  printf("  within +-2sigma of the DEUTERON rail only : %8ld  %5.1f %%\n",onlyD,100.*onlyD/nTot);
  printf("  within +-2sigma of the TRITON   rail only : %8ld  %5.1f %%\n",onlyT,100.*onlyT/nTot);
  printf("  within +-2sigma of the ALPHA    rail only : %8ld  %5.1f %%\n",onlyA,100.*onlyA/nTot);
  printf("  within +-2sigma of the PROTON   rail only : %8ld  %5.1f %%\n",onlyP,100.*onlyP/nTot);
  printf("  ambiguous (near two rails)                : %8ld  %5.1f %%\n",amb,100.*amb/nTot);
  printf("  near NO rail (background)                 : %8ld  %5.1f %%\n",none,100.*none/nTot);
  printf("\n  deuteron purity of the gated sample: %.1f %%  (of the tracks that can be assigned)\n",
         100.0*onlyD/std::max(1L,onlyD+onlyT+onlyA+onlyP));
}
