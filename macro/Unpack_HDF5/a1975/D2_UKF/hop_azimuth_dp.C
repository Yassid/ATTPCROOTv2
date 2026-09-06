/// At each hop, are the two sides on the SAME turn (a real gap -- e.g. the beam hole) or on
/// DIFFERENT turns (the ordering interleaved them)? Azimuth about the fitted circle tells them
/// apart: a gap keeps |dphi| small, an interleave gives a large |dphi|.
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>
void hop_azimuth_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                    TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=400)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   // dphi at hops, and dz at hops, for both orders
   std::vector<double> dphiWalk,dphiZ,dzWalk,dzZ;
   for(auto &j:jobs){
      if(done>=maxTracks) break;
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if(rt!=open){ if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); a=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&a); open=rt; }
      t->GetEntry(en); auto*ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
      for(auto &tr: ev->GetTrackArray()){
         if(tr.GetTrackID()!=id) continue;
         auto*hc=tr.GetHitClusterArray(); const int n=hc->size(); if(n<20) continue;
         ++done;
         double cx=tr.GetGeoCenter().first, cy=tr.GetGeoCenter().second;
         if(!std::isfinite(cx)||!std::isfinite(cy)) continue;
         std::vector<ROOT::Math::XYZPoint> st,bz;
         for(int i=0;i<n;++i) st.push_back((*hc)[i].GetPosition());
         bz=st; std::sort(bz.begin(),bz.end(),[](const ROOT::Math::XYZPoint&A,const ROOT::Math::XYZPoint&B){return A.Z()<B.Z();});
         auto scan=[&](std::vector<ROOT::Math::XYZPoint>&p, std::vector<double>&dphiOut, std::vector<double>&dzOut){
            std::vector<double> s; for(size_t i=1;i<p.size();++i) s.push_back((p[i]-p[i-1]).R());
            if(s.empty()) return;
            double mean=std::accumulate(s.begin(),s.end(),0.0)/s.size();
            for(size_t i=1;i<p.size();++i){
               if(s[i-1]<=3.0*mean) continue;
               double a1=std::atan2(p[i].Y()-cy,p[i].X()-cx), a0=std::atan2(p[i-1].Y()-cy,p[i-1].X()-cx);
               double d=a1-a0; while(d>M_PI)d-=2*M_PI; while(d<=-M_PI)d+=2*M_PI;
               dphiOut.push_back(std::fabs(d)*TMath::RadToDeg());
               dzOut.push_back(std::fabs(p[i].Z()-p[i-1].Z()));
            }
         };
         scan(st,dphiWalk,dzWalk); scan(bz,dphiZ,dzZ);
      }
   }
   if(f) f->Close();
   auto rep=[&](const char*lbl, std::vector<double>&v, std::vector<double>&z){
      if(v.empty()){printf("  %-14s none\n",lbl);return;}
      std::vector<double> s=v; std::sort(s.begin(),s.end());
      std::vector<double> zz=z; std::sort(zz.begin(),zz.end());
      int small=0,big=0; for(double d:v){ if(d<45) ++small; if(d>135) ++big; }
      printf("  %-14s n %5zu | |dphi| p10 %6.1f med %6.1f p90 %6.1f deg | <45deg %4.1f%%  >135deg %4.1f%% | |dz| med %5.2f mm\n",
             lbl, v.size(), s[s.size()/10], s[s.size()/2], s[9*s.size()/10],
             100.0*small/v.size(), 100.0*big/v.size(), zz[zz.size()/2]);
   };
   printf("\n  %d tracks\n", done);
   rep("NN walk", dphiWalk, dzWalk);
   rep("sort by z", dphiZ, dzZ);
   printf("\n  |dphi| < 45 deg  => two sides on the SAME turn: the jump is a real GAP\n");
   printf("  |dphi| > 135 deg => opposite side of the circle: turns were INTERLEAVED\n");
}
