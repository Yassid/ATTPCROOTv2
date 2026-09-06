/// Does simply SORTING CLUSTERS BY Z beat the greedy nearest-neighbour walk?
/// For a helix with axis along z, z is monotonic along the trajectory, so a z sort should give a
/// contiguous order with no hops -- and it is already what AtGenfitter does internally (:223).
/// Compares, per track: number of hops (>3x mean step) and the step distribution, both orders.
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>
static void stepStats(const std::vector<ROOT::Math::XYZPoint>&p,int&nHop,double&medStep,double&maxHop)
{
   const int n=p.size(); std::vector<double> s;
   for(int i=1;i<n;++i) s.push_back((p[i]-p[i-1]).R());
   nHop=0; maxHop=0; medStep=0; if(s.empty()) return;
   double mean=std::accumulate(s.begin(),s.end(),0.0)/s.size();
   for(double v:s) if(v>3.0*mean){++nHop; maxHop=std::max(maxHop,v);}
   std::vector<double> t=s; std::sort(t.begin(),t.end()); medStep=t[t.size()/2];
}
void zsort_vs_walk_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                      TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=800)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   long walkHopTracks=0, zHopTracks=0; long walkHops=0, zHops=0;
   std::vector<double> wMed,zMed,wMax,zMax;
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
         std::vector<ROOT::Math::XYZPoint> asStored, byZ;
         for(int i=0;i<n;++i) asStored.push_back((*hc)[i].GetPosition());
         byZ=asStored;
         std::sort(byZ.begin(),byZ.end(),[](const ROOT::Math::XYZPoint&A,const ROOT::Math::XYZPoint&B){return A.Z()<B.Z();});
         int h1,h2; double m1,m2,x1,x2;
         stepStats(asStored,h1,m1,x1); stepStats(byZ,h2,m2,x2);
         if(h1) ++walkHopTracks; if(h2) ++zHopTracks;
         walkHops+=h1; zHops+=h2; wMed.push_back(m1); zMed.push_back(m2);
         wMax.push_back(x1); zMax.push_back(x2);
      }
   }
   if(f) f->Close();
   auto med=[](std::vector<double> v){ std::sort(v.begin(),v.end()); return v.empty()?-1:v[v.size()/2]; };
   printf("\n  %d tracks\n", done);
   printf("  %-28s %12s %12s\n","","NN walk (now)","sort by z");
   printf("  %-28s %12ld %12ld\n","tracks with >=1 hop", walkHopTracks, zHopTracks);
   printf("  %-28s %12.1f%% %11.1f%%\n","  i.e.", 100.0*walkHopTracks/done, 100.0*zHopTracks/done);
   printf("  %-28s %12ld %12ld\n","total hops", walkHops, zHops);
   printf("  %-28s %12.1f %12.1f\n","median step [mm]", med(wMed), med(zMed));
   printf("  %-28s %12.1f %12.1f\n","median largest hop [mm]", med(wMax), med(zMax));
}
