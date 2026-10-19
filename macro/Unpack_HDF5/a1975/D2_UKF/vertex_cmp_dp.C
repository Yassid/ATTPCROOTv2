/// @file vertex_cmp_dp.C -- reaction-vertex distribution, data vs the g.s. sim (2026-10-12).
/// The sim beam is a pencil on axis (SetSpotRadius 0, no tilt), flat in z up to maxELoss 10 MeV. Compares, per vz
/// slice, the Spyral (fit-free) and genfit vertex x, y, r of backward protons (theta_lab 89-161, KE>1.2, data inside the
/// proton gate), and the vz distribution of all of them and of the bound-window ones.
/// Sources: gf_dp_find805 genfitter files (AtPIDEvent Spyral + AtTrackingEvent genfit), sim gs_s4001_genfitter_p.root.
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
void vertex_cmp_dp(int runMax=103){
  gSystem->Load("libAtReconstruction.so"); FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel=kFatal;
  TGraph G; { std::ifstream in("pid/proton_d2_dv1104.json"); std::string s((std::istreambuf_iterator<char>(in)),{});
    auto p=s.find("vertices"); p=s.find('[',p)+1; while(true){ auto a=s.find('[',p); if(a==std::string::npos) break; auto b=s.find(']',a); double x,y; sscanf(s.substr(a+1,b-a-1).c_str(),"%lf , %lf",&x,&y); G.AddPoint(x,y); p=b+1; if(s[s.find_first_not_of(" \n\t,",p)]==']') break; } }
  TFile fo("vertex_cmp_dp.root","RECREATE");
  TNtuple nt("v","vertices","src:run:svx:svy:svz:gvx:gvy:gvz:th:ke:ex:tcm:c2:bound");
  auto scan=[&](int src,int run,TString fn){
    TFile*f=TFile::Open(fn); if(!f||f->IsZombie()) return; TTree*t=(TTree*)f->Get("cbmsim"); TClonesArray*te=nullptr,*pe=nullptr;
    t->SetBranchAddress("AtTrackingEvent",&te); t->SetBranchAddress("AtPIDEvent",&pe);
    for(Long64_t i=0;i<t->GetEntries();++i){ if(t->GetEntry(i)<=0||!pe||pe->GetEntries()==0) continue; auto*pid=(AtPIDEvent*)pe->At(0); if(!pid) continue;
      auto*ev=(te&&te->GetEntries())?(AtTrackingEvent*)te->At(0):nullptr; std::map<int,AtFittedTrack*> fm; if(ev) for(auto&ft:ev->GetFittedTracks()) if(ft) fm[ft->GetTrackID()]=ft.get();
      for(auto&sr:pid->GetSpyral()){ if(!sr.valid) continue;
        double th=180-sr.polar*TMath::RadToDeg(), vz=1000-sr.vertex.Z(), p=299.792458*sr.brho, ke=(sqrt(p*p+938.272*938.272)-938.272)/1.023;
        if(!(th>=89&&th<=161)||!(ke>=1.2)) continue; if(src==0&&!G.IsInside(sr.sqrtdEdx,sr.brho)) continue;
        double ex=-999,tc=-999; kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc);
        double gx=-999,gy=-999,gz=-999,c2=-1; auto it=fm.find(sr.trackID);
        if(it!=fm.end()){ auto*ft=it->second; gx=ft->GetVertex().X(); gy=ft->GetVertex().Y(); gz=ft->GetVertex().Z(); double ndf=ft->GetTrackMetadata()->GetNdf(); c2=ndf>0?ft->GetTrackMetadata()->GetChi2()/ndf:1e9; }
        float row[14]={float(src),float(run),float(sr.vertex.X()),float(sr.vertex.Y()),float(vz),float(gx),float(gy),float(gz),float(th),float(ke),float(ex),float(tc),float(c2),float(ex>-1.5&&ex<0.729)};
        nt.Fill(row); } }
    f->Close(); };
  for(int r=16;r<=runMax;++r){ TString fn=Form("/mnt/f/a1975/gf_dp_find805/run_%04d_multifit_genfitter_p.root",r); if(!gSystem->AccessPathName(fn)) scan(0,r,fn); }
  scan(1,0,"/mnt/f/a1975_C16_dp_sim/gs_s4001_genfitter_p.root");
  fo.cd(); nt.Write(); fo.Close();
}
