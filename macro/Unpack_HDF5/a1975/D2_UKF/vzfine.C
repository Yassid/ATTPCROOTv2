#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
bool sel(double ke,double th,double vz,double c2,double&ex,double&tc){
  if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) return false; return kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc); }
void vzfine(){
  const int NV=22; // 50 mm bins -100..1000
  int D[3][NV]={{0}},S[3][NV]={{0}},G[3][NV]={{0}};
  auto tbin=[](double t){ return t>=10&&t<20?0: t>=20&&t<30?1: t>=30&&t<35?2:-1; };
  auto vbin=[](double v){ int b=(int)floor((v+100)/50); return (b>=0&&b<NV)?b:-1; };
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int a=tbin(tc),b=vbin(vz); if(a>=0&&b>=0) D[a][b]++;} }
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simres_dump.root"); TTree*t=(TTree*)f.Get("r"); double ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("th",&th);t->SetBranchAddress("vz",&vz);t->SetBranchAddress("c2",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int a=tbin(tc),b=vbin(vz); if(a>=0&&b>=0) S[a][b]++;} }
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simdz.root"); TTree*t=(TTree*)f.Get("s"); float vz,tc; t->SetBranchAddress("vz",&vz);t->SetBranchAddress("tc",&tc);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); int a=tbin(tc),b=vbin(vz); if(a>=0&&b>=0) G[a][b]++;} }
  const char*nm[3]={"thcm 10-20","thcm 20-30","thcm 30-35"};
  printf("vz(reco for data/sim-reco, true for gen)   per 50 mm\n  vz    ");
  for(int a=0;a<3;++a) printf("| %-10s data simR  gen ",nm[a]); printf("\n");
  for(int b=0;b<NV;++b){ printf("%5d  ",-100+50*b); for(int a=0;a<3;++a) printf("|            %4d %4d %4d ",D[a][b],S[a][b],G[a][b]); printf("\n"); }
}
