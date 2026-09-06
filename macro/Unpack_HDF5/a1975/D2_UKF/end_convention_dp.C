/// Is "highest Z = vertex end" true? Test it against the only physical marker available: the
/// vertex sits on the BEAM AXIS, so the vertex end is the one with the smaller |xy|.
/// Split by GeoTheta, because the convention should hold for forward tracks and fail for backward.
#include <cmath>
#include <fstream>
void end_convention_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                       TString gfDir="/mnt/f/a1975/gf_dp_find/", int maxTracks=800)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   // [forward/backward][cl0 nearer axis? yes/no]
   long tab[2][2]={{0,0},{0,0}};
   long dec[2][2]={{0,0},{0,0}};   // only tracks whose two ends differ by >2x in |xy|
   long ambig[2]={0,0};
   for(auto &j:jobs){
      if(done>=maxTracks) break;
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if(rt!=open){ if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); a=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&a); open=rt; }
      t->GetEntry(en); auto*ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
      for(auto &tr: ev->GetTrackArray()){
         if(tr.GetTrackID()!=id) continue;
         auto*hc=tr.GetHitClusterArray(); if(hc->size()<20) continue;
         double gt=tr.GetGeoTheta()*TMath::RadToDeg(); if(!std::isfinite(gt)) continue;
         auto c0=(*hc)[0].GetPosition(), cN=(*hc)[hc->size()-1].GetPosition();
         double xy0=std::hypot(c0.X(),c0.Y()), xyN=std::hypot(cN.X(),cN.Y());
         int bwd = (gt>90)?1:0;
         int first = (xy0<xyN)?0:1;   // 0 = cl[0] is nearer the axis (convention holds)
         tab[bwd][first]++; ++done;
         double lo=std::min(xy0,xyN), hi=std::max(xy0,xyN);
         if (lo>1e-6 && hi/lo>2.0) dec[bwd][first]++; else ambig[bwd]++;
      }
   }
   if(f) f->Close();
   printf("\n  %d tracks. Does cl[0] -- the HIGHEST-Z seed -- sit nearer the beam axis?\n\n", done);
   printf("  %-22s %14s %14s %10s\n","GeoTheta","cl[0] nearer","cl[last] nearer","convention");
   for(int b=0;b<2;++b){
      long n=tab[b][0]+tab[b][1]; if(!n) continue;
      printf("  %-22s %8ld (%3.0f%%) %8ld (%3.0f%%) %9s\n",
             b? "BACKWARD (>90 deg)":"forward (<90 deg)",
             tab[b][0],100.0*tab[b][0]/n, tab[b][1],100.0*tab[b][1]/n,
             tab[b][0]>tab[b][1]?"holds":"FAILS");
   }
   printf("\n  DECISIVE ONLY (the two ends differ by more than 2x in |xy|):\n");
   for(int b=0;b<2;++b){
      long n=dec[b][0]+dec[b][1]; if(!n) continue;
      printf("  %-22s %8ld (%3.0f%%) %8ld (%3.0f%%)   [%ld ambiguous dropped]\n",
             b? "BACKWARD (>90 deg)":"forward (<90 deg)",
             dec[b][0],100.0*dec[b][0]/n, dec[b][1],100.0*dec[b][1]/n, ambig[b]);
   }
}
