#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
/// @file fullrange_bound_dp.C -- FULL vz range with a data-driven vz efficiency f(vz) (Yassid 2026-09-24):
///   E(th,vz) = eps_sim(th,vz) x f(vz),  f = [data/eps_sim] vz profile normalised to its MAXIMUM (or a reference)
/// (selection code copied from angdist_bound_dp.C)
/// 16C(d,p)17C BOUND-state (g.s.+0.217+0.331 unresolved) angular distribution, ABSOLUTE.
/// Built 2026-09-24. Fiducial vz 200-450 mm chosen because data vz shape departs from the g.s. sim
/// beyond 450 (mid-chamber loss, cause OPEN) and 150-200 sits on the reconstruction-onset edge.
/// L_W = 316.98 mb^-1 (full D2 target, beam count x areal density, reference_a1975_luminosity)
///       x 250/971.73 (geometric; uniform gas).
/// Inputs: data  /mnt/f/a1975/caches/.explorer/exp_dp_find805.root (IC-gated slim of dp_kin_find805)
///         sim   /mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root (truth) + gs_s4001_simres_dump.root
///               (per-track reco dump from sim_res_dp.C logic, truth-closest fitted track)
///   root -l -b -q angdist_bound_dp.C
// SIM CONDITIONS both sides: keOff 0, kc off, pivot Ebeam(vz), theta_lab 89-161, KE>1.2, chi2<13.5, Ex -1.5..0.729
// acceptance eps = (sim reco-selected, reco vz in W, reco thcm bin) / (generated, TRUE vz in W, TRUE thcm bin)
const int NT=7; double TB[NT+1]={10,15,20,25,30,35,40,45};
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
bool sel(double ke,double th,double vz,double c2,double&ex,double&tc){
  if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) return false; return kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc); }
int tb(double t){ for(int i=0;i<NT;++i) if(t>=TB[i]&&t<TB[i+1]) return i; return -1; }
void fullrange_bound_dp(double Lfull=316.98,double Ltgt=971.73){
  gSystem->Load("libAtReconstruction.so");
  // generated truth
  std::vector<std::array<double,3>> gen; // tvz, tth, tke
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root"); TTree*t=(TTree*)f.Get("cbmsim"); TClonesArray*mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
    t->SetBranchStatus("*",0); t->SetBranchStatus("MCTrack*",1);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      for(int j=0;j<mc->GetEntriesFast();++j){ auto*m=(AtMCTrack*)mc->At(j); if(!m||m->GetPdgCode()!=2212) continue;
        double px=m->GetPx(),py=m->GetPy(),pz=m->GetPz(),p=sqrt(px*px+py*py+pz*pz)*1000;
        gen.push_back({m->GetStartZ()*10.0, acos(pz*1000/p)*TMath::RadToDeg(), sqrt(p*p+938.272*938.272)-938.272}); break; } } }
  printf("generated protons %zu\n",gen.size());
  // sim reco
  struct S{double ke,th,vz,c2,tke,tth,tvz;}; std::vector<S> sim;
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simres_dump.root"); TTree*t=(TTree*)f.Get("r"); S s; t->SetBranchAddress("ke",&s.ke);t->SetBranchAddress("th",&s.th);t->SetBranchAddress("vz",&s.vz);t->SetBranchAddress("c2",&s.c2);
    t->SetBranchAddress("tke",&s.tke);t->SetBranchAddress("tth",&s.tth);t->SetBranchAddress("tvz",&s.tvz);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); sim.push_back(s);} }
  // data
  struct D{double ke,th,vz,c2;}; std::vector<D> dat;
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); dat.push_back({ke,th,vz,c2});} }
  const int NZ=8; double ZB[NZ+1]={200,300,400,500,600,700,800,900,950};
  auto zb=[&](double z){ for(int k=0;k<NZ;++k) if(z>=ZB[k]&&z<ZB[k+1]) return k; return -1; };
  double Y[NT][NZ]={{0}},B[NT][NZ]={{0}},Nr[NT][NZ]={{0}},Ng[NT][NZ]={{0}};
  for(auto&d:dat){ double ex,tc; int z=zb(d.vz); if(z<0) continue; if(!sel(d.ke,d.th,d.vz,d.c2,ex,tc)) continue; int b=tb(tc); if(b<0) continue;
    if(ex>=-1.5&&ex<=0.729) Y[b][z]++; else if(ex>=-4.0&&ex<-1.5) B[b][z]+=2.229/2.5; }
  for(auto&s:sim){ double ex,tc; int z=zb(s.vz); if(z<0) continue; if(!sel(s.ke,s.th,s.vz,s.c2,ex,tc)) continue; if(ex<-1.5||ex>0.729) continue; int b=tb(tc); if(b>=0) Nr[b][z]++; }
  for(auto&g:gen){ int z=zb(g[0]); if(z<0) continue; double ex,tc; if(!kine2b(Eb(g[0]),g[1]*TMath::DegToRad(),g[2],ex,tc)) continue; int b=tb(tc); if(b>=0) Ng[b][z]++; }
  auto eps=[&](int b,int z){ return Ng[b][z]>0?Nr[b][z]/Ng[b][z]:0.; };
  printf("\nMAP: data counts (net) / eps_sim / eps_sim cell, then data/eps_sim PER 100 mm -- rows thcm, cols vz\n  thcm   ");
  for(int z=0;z<NZ;++z) printf("   %3.0f-%3.0f      ",ZB[z],ZB[z+1]); printf("\n");
  for(int b=1;b<=5;++b){ printf("  %2.0f-%2.0f ",TB[b],TB[b+1]); for(int z=0;z<NZ;++z){ double e=eps(b,z); printf("  %3.0f/%.2f=%5.1f ",Y[b][z]-B[b][z],e,e>0?(Y[b][z]-B[b][z])/e*100/(ZB[z+1]-ZB[z]):0);} printf("\n"); }
  // corrected density per mm in each cell (recommended theta bins 1..5 only for f)
  double Dz[3][NZ]={{0}}; // 0 all thcm 15-40, 1 thcm 15-35, 2 thcm 35-40
  for(int z=0;z<NZ;++z) for(int b=1;b<=5;++b){ double e=eps(b,z); if(e<=0) continue; double v=(Y[b][z]-B[b][z])/e/(ZB[z+1]-ZB[z]);
    Dz[0][z]+=v; if(b<=4) Dz[1][z]+=v; else Dz[2][z]+=v; }
  printf("\nDATA/eps_sim PER MM vs vz (the data-only vz efficiency before normalisation)\n  vz        all15-40   15-35    35-40   | eps_sim(thcm 25-30)  eps_sim(35-40)\n");
  for(int z=0;z<NZ;++z) printf("  %3.0f-%3.0f   %6.3f   %6.3f   %6.3f   |   %.2f               %.2f\n",ZB[z],ZB[z+1],Dz[0][z],Dz[1][z],Dz[2][z],eps(3,z),eps(5,z));
  auto norm=[&](int g,int mode){ std::vector<double> f(NZ); double ref=0;
    if(mode==0){ for(int z=0;z<NZ;++z) ref=std::max(ref,Dz[g][z]); }                       // MAXIMUM
    else if(mode==1){ for(int z=0;z<5;++z) ref=std::max(ref,Dz[g][z]); }                    // max over 200-700 only
    else { ref=(Dz[g][0]+Dz[g][1]+0.5*Dz[g][2])/2.5; }                                     // 200-450 plateau
    for(int z=0;z<NZ;++z) f[z]=ref>0?Dz[g][z]/ref:0; return f; };
  const char*MN[3]={"max (all vz)","max over 200-700","200-450 plateau"};
  double Lr=Lfull*(ZB[NZ]-ZB[0])/Ltgt;
  printf("\nL for 150-940 = %.2f mb^-1\n",Lr);
  // results: rows theta, columns method
  struct Res{double v[NT];};
  std::vector<std::pair<std::string,Res>> out;
  for(int mode=0;mode<3;++mode) for(int grp=0;grp<2;++grp){ // grp0: one f from all theta; grp1: f per theta group
    auto fA=norm(0,mode), f1=norm(1,mode), f2=norm(2,mode); Res r;
    if(mode==0&&grp==0){ printf("\nf(vz) normalised to %s:\n  vz      ",MN[mode]); for(int z=0;z<NZ;++z) printf(" %4.0f",ZB[z]); printf("\n  all    "); for(int z=0;z<NZ;++z) printf(" %4.2f",fA[z]); printf("\n  15-35  "); for(int z=0;z<NZ;++z) printf(" %4.2f",f1[z]); printf("\n  35-40  "); for(int z=0;z<NZ;++z) printf(" %4.2f",f2[z]); printf("\n"); }
    for(int b=0;b<NT;++b){ double sum=0; double dO=2*TMath::Pi()*(cos(TB[b]*TMath::DegToRad())-cos(TB[b+1]*TMath::DegToRad()));
      for(int z=0;z<NZ;++z){ double e=eps(b,z); double f= grp==0? fA[z] : (b<=4? f1[z]: f2[z]); if(e<=0||f<=0.05) continue; sum+=(Y[b][z]-B[b][z])/(e*f); }
      r.v[b]=sum/(dO*Lr); }
    out.push_back({std::string(MN[mode])+(grp?", f per theta group":", one f"),r}); }
  double W[NT]={0.148,1.799,3.023,5.232,6.707,3.303,0.364};
  printf("\nFULL-RANGE dsigma/dOmega [mb/sr] vs the vz 200-450 window result\n  thcm   window  ");
  for(auto&o:out) printf(" | %-28s",o.first.c_str()); printf("\n");
  for(int b=1;b<=5;++b){ printf("  %2.0f-%2.0f  %5.2f  ",TB[b],TB[b+1],W[b]); for(auto&o:out) printf(" |  %5.2f  (x%.2f)                ",o.second.v[b],o.second.v[b]/W[b]); printf("\n"); }
}
