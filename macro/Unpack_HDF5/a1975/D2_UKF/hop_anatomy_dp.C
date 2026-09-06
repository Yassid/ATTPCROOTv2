/// WHERE do the greedy-walk hops happen and WHAT is on the far side of them?
/// If hops sit at the end of the array and the piece after is small, the walk is not
/// mis-ordering the track -- it is APPENDING clusters that are not on it.
#include <cmath>
#include <vector>
#include <fstream>
void hop_anatomy_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                    TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=400, double jumpFac=3.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   std::vector<double> firstHopFrac, tailFrac, tailQFrac, tailXY, mainXY;
   long nTr=0, nWithHop=0;
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
         ++done; ++nTr;
         std::vector<double> steps(n,0); double mean=0;
         for(int i=1;i<n;++i){ steps[i]=((*hc)[i].GetPosition()-(*hc)[i-1].GetPosition()).R(); mean+=steps[i]; }
         mean/=(n-1);
         int first=-1;
         for(int i=1;i<n;++i) if(steps[i]>jumpFac*mean){ first=i; break; }
         if(first<0) continue;
         ++nWithHop;
         firstHopFrac.push_back((double)first/n);
         tailFrac.push_back((double)(n-first)/n);
         double qTail=0,qAll=0,xyTail=0,xyMain=0;
         for(int i=0;i<n;++i){ double q=(*hc)[i].GetCharge(); qAll+=q;
            auto p=(*hc)[i].GetPosition(); double xy=std::hypot(p.X(),p.Y());
            if(i>=first){ qTail+=q; xyTail+=xy; } else xyMain+=xy; }
         tailQFrac.push_back(qAll>0?qTail/qAll:0);
         if(n-first>0) tailXY.push_back(xyTail/(n-first));
         if(first>0)   mainXY.push_back(xyMain/first);
      }
   }
   if(f) f->Close();
   auto med=[](std::vector<double> v){ if(v.empty()) return -1.0; std::sort(v.begin(),v.end()); return v[v.size()/2]; };
   auto qq=[](std::vector<double> v,double p){ if(v.empty()) return -1.0; std::sort(v.begin(),v.end()); return v[(size_t)(p*(v.size()-1))]; };
   printf("\n  %ld tracks, %ld with a hop\n", nTr, nWithHop);
   printf("  position of the FIRST hop (fraction of array): p10 %.2f  median %.2f  p90 %.2f\n",
          qq(firstHopFrac,.1), med(firstHopFrac), qq(firstHopFrac,.9));
   printf("  size of the piece AFTER it   (fraction)      : p10 %.2f  median %.2f  p90 %.2f\n",
          qq(tailFrac,.1), med(tailFrac), qq(tailFrac,.9));
   printf("  CHARGE fraction in that piece                : p10 %.3f  median %.3f  p90 %.3f\n",
          qq(tailQFrac,.1), med(tailQFrac), qq(tailQFrac,.9));
   printf("  mean |xy| before the hop %.1f mm   after it %.1f mm\n", med(mainXY), med(tailXY));
}
