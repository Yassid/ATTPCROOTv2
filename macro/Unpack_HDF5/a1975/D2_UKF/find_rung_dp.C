/// Which rung does FIND actually land on, and does it depend on track length?
/// The stored fit's number of smoothed positions over the track's cluster count recovers the
/// fraction FIND kept -- the only trace it leaves, since findPct is not written to the output.
/// Tests Yassid's proposal: if long tracks always land low, the high rungs are wasted on them and
/// the range could be made length-dependent (lower floor at the same cost).
#include <cmath>
#include <vector>
#include <map>
void find_rung_dp(TString runsCSV="0016,0020,0026,0034",
                  TString gfDir="/mnt/f/a1975/gf_dp_findscan/", int minCl=20)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TObjArray *toks=runsCSV.Tokenize(",");
   // frac kept, bucketed by cluster count
   const int NB=4; const int lo[NB]={20,100,300,600}, hi[NB]={100,300,600,100000};
   std::vector<double> frac[NB]; long full[NB]={0,0,0,0}, tot[NB]={0,0,0,0};
   for(int r=0;r<toks->GetEntries();++r){
      TString num=((TObjString*)toks->At(r))->GetString();
      auto*f=TFile::Open(gfDir+"run_"+num+"_multifit_genfitter_p.root");
      if(!f||f->IsZombie()) continue;
      auto*t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
      t->SetBranchAddress("AtTrackingEvent",&a);
      for(Long64_t i=0;i<t->GetEntries();++i){
         t->GetEntry(i); auto*ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
         std::map<int,int> ncl;
         for(auto &tr: ev->GetTrackArray()) ncl[tr.GetTrackID()]=tr.GetHitClusterArray()->size();
         for(auto &ft: ev->GetFittedTracks()){
            auto it=ncl.find(ft->GetTrackID()); if(it==ncl.end()) continue;
            const int n=it->second; if(n<minCl) continue;
            const int pts=ft->GetSmoothedPositions().size();
            if(pts<3) continue;
            double fr=(double)pts/n;
            int b=-1; for(int k=0;k<NB;++k) if(n>=lo[k]&&n<hi[k]) b=k;
            if(b<0) continue;
            ++tot[b]; frac[b].push_back(fr);
            if(fr>0.93) ++full[b];      // effectively full length
         }
      }
      f->Close();
   }
   printf("\n  fraction of the track FIND kept (stored fit points / clusters)\n\n");
   printf("  %-14s %8s %10s %9s %9s %9s %9s\n","clusters","n","kept full","p10","median","p90","min");
   for(int b=0;b<NB;++b){
      if(frac[b].empty()) continue;
      std::sort(frac[b].begin(),frac[b].end());
      auto q=[&](double p){ return frac[b][(size_t)(p*(frac[b].size()-1))]; };
      printf("  %5d-%-8d %8ld %9.0f%% %9.2f %9.2f %9.2f %9.2f\n",
             lo[b], hi[b]>99999?9999:hi[b], tot[b], 100.0*full[b]/tot[b], q(0.1), q(0.5), q(0.9), frac[b][0]);
   }
   printf("\n  'kept full' = fraction of tracks FIND left at (nearly) 100 %%.\n");
   printf("  If long tracks rarely stay full and land low, the high rungs are wasted on them.\n");
}
