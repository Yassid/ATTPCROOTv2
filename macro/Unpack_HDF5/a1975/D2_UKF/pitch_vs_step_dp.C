/// Why the 3D nearest-neighbour cluster walk strands clusters: the z advance per TURN (the helix
/// pitch) compared with the cluster-to-cluster step. When pitch/step is small the nearest cluster
/// in 3D lies on the NEXT TURN, so the walk cuts across and comes back later (a "hop").
#include <cmath>
#include <vector>
#include <fstream>
void pitch_vs_step_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                      TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=800,
                      double jumpFac=3.0,
                      TString csvOut="/tmp/claude-1000/-home-yassid/d509a965-524d-4548-9401-35c7b60d6646/scratchpad/pitch.csv")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   std::ofstream o(csvOut.Data()); o<<"run,entry,tid,ncl,R,meanStep,pitch,pitchOverStep,nHop,maxHop\n";
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
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
         double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
         for(int i=0;i<n;++i){auto p=(*hc)[i].GetPosition();double x=p.X(),y=p.Y();
            Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y;}
         double C=n*Sxx-Sx*Sx,D=n*Sxy-Sx*Sy,E=n*Sxxx+n*Sxyy-(Sxx+Syy)*Sx;
         double G=n*Syy-Sy*Sy,H=n*Sxxy+n*Syyy-(Sxx+Syy)*Sy,den=2*(C*G-D*D);
         if(std::fabs(den)<1e-9) continue;
         double cx=(E*G-D*H)/den, cy=(C*H-D*E)/den, R=0;
         for(int i=0;i<n;++i){auto p=(*hc)[i].GetPosition(); R+=std::hypot(p.X()-cx,p.Y()-cy);} R/=n;
         // per-step |dz| and |dphi| along the STORED order, medians -> pitch = 2pi * median(dz/dphi)
         std::vector<double> steps, ratio;
         for(int i=1;i<n;++i){
            auto p=(*hc)[i].GetPosition(), q=(*hc)[i-1].GetPosition();
            steps.push_back((p-q).R());
            double a1=std::atan2(p.Y()-cy,p.X()-cx), a0=std::atan2(q.Y()-cy,q.X()-cx);
            double dphi=a1-a0; while(dphi>M_PI)dphi-=2*M_PI; while(dphi<=-M_PI)dphi+=2*M_PI;
            if(std::fabs(dphi)>1e-3) ratio.push_back(std::fabs((p.Z()-q.Z())/dphi));
         }
         if(steps.empty()||ratio.empty()) continue;
         std::vector<double> s2=steps; std::sort(s2.begin(),s2.end());
         double medStep=s2[s2.size()/2];
         std::sort(ratio.begin(),ratio.end());
         double pitch=2*M_PI*ratio[ratio.size()/2];
         double meanStep=0; for(double v:steps) meanStep+=v; meanStep/=steps.size();
         int nHop=0; double maxHop=0;
         for(double v:steps){ if(v>jumpFac*meanStep){++nHop; maxHop=std::max(maxHop,v);} }
         o<<rt<<","<<en<<","<<id<<","<<n<<","<<R<<","<<medStep<<","<<pitch<<","
          <<(medStep>0?pitch/medStep:-1)<<","<<nHop<<","<<maxHop<<"\n";
         ++done;
      }
   }
   if(f) f->Close(); o.close();
   printf("  %d tracks -> %s\n", done, csvOut.Data());
}
