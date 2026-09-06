/// @file zuturn_dp.C
/// @brief How often does the stored cluster array U-TURN in z, and on what kind of track?
///
/// A greedy walk must consume every cluster, so any it skips get APPENDED AT THE END -- the tail
/// of the array becomes "everything the walk missed" rather than the far end of the track. The
/// signature is z going down and then back up (run_0016 e31316: 683 -> 122 -> 649).
///
/// METRIC: total variation of z along the array, divided by the z span.
///   1.0 = monotonic traversal (what a correct ordering gives)
///   2.0 = one full U-turn (down and all the way back)
/// Correlated against turns swept and cluster count, to answer why only some events are affected.
#include <cmath>
#include <vector>
#include <fstream>
void zuturn_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
               TString gfDir="/mnt/f/a1975/gf_dp_findscan/", int maxTracks=800)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   struct Rec { double tv; double turns; int ncl; };
   std::vector<Rec> v;
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
         double tv=0, zlo=1e9, zhi=-1e9;
         for(int i=0;i<n;++i){ double z=(*hc)[i].GetPosition().Z(); zlo=std::min(zlo,z); zhi=std::max(zhi,z);
            if(i) tv += std::fabs(z-(*hc)[i-1].GetPosition().Z()); }
         if(zhi-zlo < 20) continue;
         double cx=tr.GetGeoCenter().first, cy=tr.GetGeoCenter().second;
         double cum=0, prev=0; bool first=true;
         if(std::isfinite(cx)&&std::isfinite(cy))
            for(int i=0;i<n;++i){ auto p=(*hc)[i].GetPosition();
               double ph=std::atan2(p.Y()-cy,p.X()-cx);
               if(!first){ double d=ph-prev; while(d>M_PI)d-=2*M_PI; while(d<=-M_PI)d+=2*M_PI; cum+=d; }
               prev=ph; first=false; }
         v.push_back({tv/(zhi-zlo), std::fabs(cum)/(2*M_PI), n});
         ++done;
      }
   }
   if(f) f->Close();
   std::vector<double> tvs; for(auto&x:v) tvs.push_back(x.tv);
   std::sort(tvs.begin(),tvs.end());
   long bad=0; for(auto&x:v) if(x.tv>1.5) ++bad;
   printf("\n  %zu tracks.  z total-variation / z span   (1.0 = monotonic, 2.0 = one U-turn)\n", v.size());
   printf("    p10 %.2f   median %.2f   p75 %.2f   p90 %.2f   max %.2f\n",
          tvs[tvs.size()/10], tvs[tvs.size()/2], tvs[3*tvs.size()/4], tvs[9*tvs.size()/10], tvs.back());
   printf("    tracks with ratio > 1.5 (a real U-turn): %ld = %.1f%%\n", bad, 100.0*bad/v.size());
   printf("\n  by TURNS SWEPT:\n    %-14s %6s %10s %12s\n","turns","n","median tv","frac >1.5");
   const double tlo[5]={0,1,2,4,8}, thi[5]={1,2,4,8,1e9};
   for(int b=0;b<5;++b){ std::vector<double> s; long nb=0;
      for(auto&x:v) if(x.turns>=tlo[b]&&x.turns<thi[b]){ s.push_back(x.tv); if(x.tv>1.5) ++nb; }
      if(s.size()<5) continue; std::sort(s.begin(),s.end());
      printf("    %5.0f - %-6.0f %6zu %10.2f %11.0f%%\n", tlo[b], thi[b]>1e8?99:thi[b], s.size(),
             s[s.size()/2], 100.0*nb/s.size()); }
   printf("\n  by CLUSTER COUNT:\n    %-14s %6s %10s %12s\n","clusters","n","median tv","frac >1.5");
   const int clo[4]={20,100,300,600}, chi[4]={100,300,600,100000};
   for(int b=0;b<4;++b){ std::vector<double> s; long nb=0;
      for(auto&x:v) if(x.ncl>=clo[b]&&x.ncl<chi[b]){ s.push_back(x.tv); if(x.tv>1.5) ++nb; }
      if(s.size()<5) continue; std::sort(s.begin(),s.end());
      printf("    %5d - %-6d %6zu %10.2f %11.0f%%\n", clo[b], chi[b]>99999?9999:chi[b], s.size(),
             s[s.size()/2], 100.0*nb/s.size()); }
}
