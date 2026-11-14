#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
/// @file angdist_bound_dp_pertheta.C -- PER-THETA vz fiducial ranges (Yassid 2026-09-24): each theta_cm bin uses
/// the vz range where data/eps_sim is flat (fullrange_bound_dp.C map), with its own eps and its own L length.
/// (selection code copied from angdist_bound_dp.C; header below describes the single-window original)
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
void angdist_bound_dp_pertheta(double Lfull=316.98,double Ltgt=971.73,const char*csv="dp_bound_angdist_pertheta.csv"){
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
  //                   10-15 15-20 20-25 25-30 30-35 35-40 40-45
  double ZL[NT]={        200,  200,  200,  200,  200,  200,  200};
  double ZH[NT]={        400,  500,  400,  400,  400,  950,  400};
  int    USE[NT]={         0,    1,    1,    1,    1,    1,    0};
  // flatness evidence: each 100/50 mm cell vs the bin's 200-400 reference, in data/eps_sim per 100 mm
  const int NZ=8; double ZB[NZ+1]={200,300,400,500,600,700,800,900,950};
  auto zb=[&](double z){ for(int k=0;k<NZ;++k) if(z>=ZB[k]&&z<ZB[k+1]) return k; return -1; };
  double Yc[NT][NZ]={{0}},Bc[NT][NZ]={{0}},Rc[NT][NZ]={{0}},Gc[NT][NZ]={{0}};
  for(auto&d:dat){ double ex,tc; int z=zb(d.vz); if(z<0) continue; if(!sel(d.ke,d.th,d.vz,d.c2,ex,tc)) continue; int b=tb(tc); if(b<0) continue;
    if(ex>=-1.5&&ex<=0.729) Yc[b][z]++; else if(ex>=-4.0&&ex<-1.5) Bc[b][z]+=2.229/2.5; }
  for(auto&s:sim){ double ex,tc; int z=zb(s.vz); if(z<0) continue; if(!sel(s.ke,s.th,s.vz,s.c2,ex,tc)) continue; if(ex<-1.5||ex>0.729) continue; int b=tb(tc); if(b>=0) Rc[b][z]++; }
  for(auto&g:gen){ int z=zb(g[0]); if(z<0) continue; double ex,tc; if(!kine2b(Eb(g[0]),g[1]*TMath::DegToRad(),g[2],ex,tc)) continue; int b=tb(tc); if(b>=0) Gc[b][z]++; }
  printf("\nFLATNESS: pull of each cell vs the bin's 200-400 reference (data/eps_sim per mm); [x] = inside the chosen range\n  thcm   ref/100mm ");
  for(int z=0;z<NZ;++z) printf("  %3.0f-%3.0f ",ZB[z],ZB[z+1]); printf("\n");
  for(int b=0;b<NT;++b){ double y=Yc[b][0]+Yc[b][1]-Bc[b][0]-Bc[b][1], e=(Gc[b][0]+Gc[b][1])>0?(Rc[b][0]+Rc[b][1])/(Gc[b][0]+Gc[b][1]):0; double ref=e>0?y/e/200.:0;
    printf("  %2.0f-%2.0f  %6.1f   ",TB[b],TB[b+1],ref*100);
    for(int z=0;z<NZ;++z){ double ez=Gc[b][z]>0?Rc[b][z]/Gc[b][z]:0, len=ZB[z+1]-ZB[z], exp=ref*len*ez, obs=Yc[b][z]-Bc[b][z];
      double pull= exp>0? (obs-exp)/sqrt(std::max(1.,exp)) : 0; bool in=ZB[z]>=ZL[b]&&ZB[z+1]<=ZH[b];
      printf("  %s%+5.1f%s ",in?"[":" ",pull,in?"]":" "); }
    printf("\n"); }
  // the result
  FILE*fo=fopen(csv,"w");
  fprintf(fo,"# 16C(d,p)17C bound states (g.s.+0.217+0.331), PER-THETA vz fiducial ranges; L_b = %.2f x (vz length)/%.2f mb^-1; cm frame\n",Lfull,Ltgt);
  fprintf(fo,"# cuts: theta_lab 89-161, KE>1.2, chi2/ndf<13.5, Ex -1.5..0.729 (pivot Ebeam(vz)=184.25-0.01103(vz-500)), keOff 0, kc off; bg = sideband Ex -4..-1.5 scaled 2.229/2.5\n");
  fprintf(fo,"# eps = sim reco-selected (reco vz in range, reco thcm bin) / generated (true vz in range, true thcm bin). errors: stat (data + bg + eps binomial), NO systematics\n");
  fprintf(fo,"thcm_lo,thcm_hi,thcm_mid,vz_lo,vz_hi,L_mb-1,Y,bg,eps,eps_err,dsdo_mb_sr,err_mb_sr,use\n");
  double W[NT]={0.148,1.799,3.023,5.232,6.707,3.303,0.364};
  printf("\n  thcm   vz range   L[mb^-1]    Y     bg    eps(+-)          dsigma/dOmega [mb/sr]   | vz200-450 window\n");
  for(int b=0;b<NT;++b){ double z0=ZL[b],z1=ZH[b],LW=Lfull*(z1-z0)/Ltgt, Y=0,Bg=0,Nr=0,Ng=0;
    for(auto&d:dat){ double ex,tc; if(d.vz<z0||d.vz>=z1) continue; if(!sel(d.ke,d.th,d.vz,d.c2,ex,tc)) continue; if(tb(tc)!=b) continue;
      if(ex>=-1.5&&ex<=0.729) Y++; else if(ex>=-4.0&&ex<-1.5) Bg++; }
    for(auto&s:sim){ double ex,tc; if(s.vz<z0||s.vz>=z1) continue; if(!sel(s.ke,s.th,s.vz,s.c2,ex,tc)) continue; if(ex<-1.5||ex>0.729) continue; if(tb(tc)==b) Nr++; }
    for(auto&g:gen){ if(g[0]<z0||g[0]>=z1) continue; double ex,tc; if(!kine2b(Eb(g[0]),g[1]*TMath::DegToRad(),g[2],ex,tc)) continue; if(tb(tc)==b) Ng++; }
    double e=Ng>0?Nr/Ng:0, ee=Ng>0?sqrt(e*(1-e)/Ng):0, dO=2*TMath::Pi()*(cos(TB[b]*TMath::DegToRad())-cos(TB[b+1]*TMath::DegToRad()));
    double bgs=Bg*(2.229/2.5), net=Y-bgs, v=(e>0)?net/(e*dO*LW):0;
    double rel2=(net>0)?(Y+bgs*bgs/std::max(1.,Bg))/(net*net):0; rel2+=e>0?(ee/e)*(ee/e):0; double er=v*sqrt(rel2);
    printf("  %2.0f-%2.0f  %3.0f-%3.0f   %6.2f   %4.0f  %5.1f   %.3f(%.3f)    %6.3f +- %.3f %s  | %.3f\n",TB[b],TB[b+1],z0,z1,LW,Y,bgs,e,ee,v,er,USE[b]?"":"(not rec.)",W[b]);
    fprintf(fo,"%.0f,%.0f,%.1f,%.0f,%.0f,%.3f,%.0f,%.2f,%.4f,%.4f,%.4f,%.4f,%d\n",TB[b],TB[b+1],0.5*(TB[b]+TB[b+1]),z0,z1,LW,Y,bgs,e,ee,v,er,USE[b]); }
  fclose(fo); printf("wrote %s\n",csv);
}
