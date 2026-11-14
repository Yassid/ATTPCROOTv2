#include <TClonesArray.h>
#include <TSystem.h>
#include "AtMCTrack.h"
#include "AtMCPoint.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
// per generated g.s. proton: truth vz, thcm, and z extent of its MC points (mm, sim frame)
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
void simdz(){
  gSystem->Load("libAtSimulationData.so");
  TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root"); TTree*t=(TTree*)f.Get("cbmsim");
  TClonesArray*mc=nullptr,*pt=nullptr; t->SetBranchAddress("MCTrack",&mc); t->SetBranchAddress("AtTpcPoint",&pt);
  TFile fo("/mnt/f/a1975_C16_dp_sim/gs_s4001_simdz.root","RECREATE"); TTree o("s","s"); float vz,tc,ke,thl,zmin,zmax,zlast,rmax; int np;
  o.Branch("vz",&vz);o.Branch("tc",&tc);o.Branch("ke",&ke);o.Branch("thl",&thl);o.Branch("zmin",&zmin);o.Branch("zmax",&zmax);o.Branch("zlast",&zlast);o.Branch("rmax",&rmax);o.Branch("np",&np);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); int jp=-1;
    for(int j=0;j<mc->GetEntriesFast();++j){auto*m=(AtMCTrack*)mc->At(j); if(m&&m->GetPdgCode()==2212){jp=j;break;}}
    if(jp<0) continue; auto*m=(AtMCTrack*)mc->At(jp);
    double px=m->GetPx(),py=m->GetPy(),pz=m->GetPz(),p=sqrt(px*px+py*py+pz*pz)*1000;
    vz=m->GetStartZ()*10; thl=acos(pz*1000/p)*TMath::RadToDeg(); ke=sqrt(p*p+938.272*938.272)-938.272;
    double ex,t2; if(!kine2b(Eb(vz),thl*TMath::DegToRad(),ke,ex,t2)) continue; tc=t2;
    zmin=1e9;zmax=-1e9;zlast=0;rmax=0;np=0; double tl=-1;
    for(int k=0;k<pt->GetEntriesFast();++k){auto*q=(AtMCPoint*)pt->At(k); if(q->GetTrackID()!=jp) continue; np++;
      double z=q->GetZ()*10; zmin=std::min<double>(zmin,z); zmax=std::max<double>(zmax,z);
      double r=10*hypot(q->GetX(),q->GetY()); rmax=std::max<double>(rmax,r); if(q->GetTime()>tl){tl=q->GetTime(); zlast=z;} }
    o.Fill(); }
  o.Write(); printf("filled %lld\n",o.GetEntries());
  // table: thcm bins x vz slices, median (zlast - vz)
  double TB[]={10,15,20,25,30,35,40,45}; double VZ[]={0,200,400,600,800,950};
  printf("median zlast-vz [mm] (sim frame) / n\n thcm  "); for(int v=0;v<5;++v) printf("  vz%3.0f-%3.0f ",VZ[v],VZ[v+1]); printf("\n");
  for(int b=0;b<7;++b){ printf("%2.0f-%2.0f ",TB[b],TB[b+1]);
    for(int v=0;v<5;++v){ std::vector<double> d; for(Long64_t e=0;e<o.GetEntries();++e){o.GetEntry(e); if(tc>=TB[b]&&tc<TB[b+1]&&vz>=VZ[v]&&vz<VZ[v+1]&&np>0) d.push_back(zlast-vz);}
      if(d.size()){std::sort(d.begin(),d.end()); printf(" %7.0f/%-4zu ",d[d.size()/2],d.size());} else printf("      -     "); } printf("\n"); }
}
