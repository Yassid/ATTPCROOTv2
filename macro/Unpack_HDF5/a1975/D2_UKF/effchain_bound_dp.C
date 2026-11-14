/// @file effchain_bound_dp.C -- ABSOLUTE efficiency chain for 16C(d,p) bound-state protons, data vs g.s. sim.
/// Spyral (fit-free) kinematics exist for EVERY pattern track; genfit only for gated (data) / all (sim).
/// Mapping measured on run_0016 (2818 tracks): theta_fit = 180 - polar_spy (+1.7 deg), vz_fit = 1000 - z_spy
/// (+2.3 mm), KE_spy/KE_fit = 1.023.
///   N0 = valid Spyral track whose SPYRAL kinematics give Ex in the bound window (theta_cm bin, vz range)
///   N1 = N0 inside the proton gate pid/proton_d2_dv1104.json
///   N2 = N1 with a genfit fit that passes the analysis cuts AND fitted Ex in the bound window
/// No IC gate (not in the genfitter files).  Built 2026-09-24.
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
const int NT=7; double TB[NT+1]={10,15,20,25,30,35,40,45};
double ZL[NT]={200,200,200,200,200,200,200}, ZH[NT]={400,500,400,400,400,950,400};
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
int tb(double t){ for(int i=0;i<NT;++i) if(t>=TB[i]&&t<TB[i+1]) return i; return -1; }
TGraph* gate(){ auto g=new TGraph; std::ifstream in("pid/proton_d2_dv1104.json"); std::string s((std::istreambuf_iterator<char>(in)),{});
  auto p=s.find("vertices"); p=s.find('[',p)+1; while(true){ auto a=s.find('[',p); if(a==std::string::npos) break; auto b=s.find(']',a); double x,y; sscanf(s.substr(a+1,b-a-1).c_str(),"%lf , %lf",&x,&y); g->AddPoint(x,y); p=b+1; if(s[s.find_first_not_of(" \n\t,",p)]==']') break; } return g; }
double N[2][3][NT]; TNtuple* pts[2];
void scan(int src,TString fn,TGraph*G){
  TFile*f=TFile::Open(fn); if(!f||f->IsZombie()) return; TTree*t=(TTree*)f->Get("cbmsim"); TClonesArray*te=nullptr,*pe=nullptr;
  t->SetBranchAddress("AtTrackingEvent",&te); t->SetBranchAddress("AtPIDEvent",&pe);
  for(Long64_t i=0;i<t->GetEntries();++i){ if(t->GetEntry(i)<=0) continue; if(!pe||pe->GetEntries()==0) continue;
    auto*pid=(AtPIDEvent*)pe->At(0); if(!pid) continue; auto*ev=(te&&te->GetEntries())?(AtTrackingEvent*)te->At(0):nullptr;
    std::map<int,AtFittedTrack*> fm; if(ev) for(auto&ft:ev->GetFittedTracks()) if(ft) fm[ft->GetTrackID()]=ft.get();
    for(auto&sr:pid->GetSpyral()){ if(!sr.valid) continue;
      double th=180-sr.polar*TMath::RadToDeg(), vz=1000-sr.vertex.Z(), p=299.792458*sr.brho, ke=(sqrt(p*p+938.272*938.272)-938.272)/1.023;
      if(!(th>=89&&th<=161)||!(ke>=1.2)) continue;
      double ex,tc; if(!kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc)) continue; if(ex<-1.5||ex>0.729) continue;
      int b=tb(tc); if(b<0||vz<ZL[b]||vz>=ZH[b]) continue;
      bool in=G->IsInside(sr.sqrtdEdx,sr.brho);
      N[src][0][b]++; if(in) N[src][1][b]++;
      pts[src]->Fill(b,sr.sqrtdEdx,sr.brho,in,vz);
      if(!in) continue;
      auto it=fm.find(sr.trackID); if(it==fm.end()) continue; auto*ft=it->second; auto&k=ft->GetKinematicsXtr();
      double kf=k.kineticEnergy, tf=k.theta*TMath::RadToDeg(), vf=ft->GetVertex().Z(), ndf=ft->GetTrackMetadata()->GetNdf(), c2=ndf>0?ft->GetTrackMetadata()->GetChi2()/ndf:1e9;
      if(!(c2<=13.5)||!(tf>=89&&tf<=161)||!(kf>=1.2)) continue; double exf,tcf; if(!kine2b(Eb(vf),tf*TMath::DegToRad(),kf,exf,tcf)) continue;
      if(exf<-1.5||exf>0.729) continue; N[src][2][b]++; } }
  f->Close(); }
void effchain_bound_dp(){
  gSystem->Load("libAtReconstruction.so"); FairLogger::GetLogger()->SetLogScreenLevel("FATAL"); gErrorIgnoreLevel=kFatal;
  auto G=gate(); printf("gate vertices %d\n",G->GetN());
  pts[0]=new TNtuple("pd","data N0","b:sq:brho:in:vz"); pts[1]=new TNtuple("ps","sim N0","b:sq:brho:in:vz");
  for(int r=16;r<=103;++r){ TString fn=Form("/mnt/f/a1975/gf_dp_find805/run_%04d_multifit_genfitter_p.root",r); if(gSystem->AccessPathName(fn)) continue; scan(0,fn,G); }
  scan(1,"/mnt/f/a1975_C16_dp_sim/gs_s4001_genfitter_p.root",G);
  printf("\n  thcm   vz       | DATA  N0    gate N1/N0    fit+cuts N2/N1   N2/N0  | SIM  N0    gate N1/N0    fit+cuts N2/N1   N2/N0 | data/sim N2/N0\n");
  for(int b=0;b<NT;++b){ double r[2][3];
    for(int s=0;s<2;++s){ r[s][0]=N[s][0][b]; r[s][1]=N[s][0][b]>0?N[s][1][b]/N[s][0][b]:0; r[s][2]=N[s][1][b]>0?N[s][2][b]/N[s][1][b]:0; }
    printf("  %2.0f-%2.0f %3.0f-%3.0f | %5.0f   %.3f          %.3f            %.3f  | %5.0f   %.3f          %.3f            %.3f | %.2f\n",TB[b],TB[b+1],ZL[b],ZH[b],
      r[0][0],r[0][1],r[0][2],r[0][1]*r[0][2],r[1][0],r[1][1],r[1][2],r[1][1]*r[1][2],(r[1][1]*r[1][2]>0)?(r[0][1]*r[0][2])/(r[1][1]*r[1][2]):0); }
  TFile fo("effchain_bound_dp.root","RECREATE"); pts[0]->Write(); pts[1]->Write(); G->Write("gate"); fo.Close();
}
