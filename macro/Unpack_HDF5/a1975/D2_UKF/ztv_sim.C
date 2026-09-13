/// THE NOISE FLOOR of the z total-variation metric, measured where the ordering is known good.
///
/// In SIMULATION the hit array arrives already in track order (the digitiser produces hits pad by
/// pad along the trajectory; measured hit-path ratio 1.01 against 1.63 in data), so the
/// clusteriser and the greedy walk have nothing to get wrong.  Whatever tv/span the sim clusters
/// give is therefore what a CORRECT order costs given real z digitisation -- which is the number
/// the data measurement has to be read against, and which no data-only comparison can supply.
///
/// Same finder, same field, same gas as the a1975 (d,p) protons: 14C(d,p) at 2.85 T.
#include <cmath>
#include <vector>
void ztv_sim(const char *file = "/mnt/f/a1954_C14dp_sm/gs/reco.root", int maxTrk = 800)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TFile *f = TFile::Open(file);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);
   std::vector<double> tv; std::vector<int> ncl;
   int done = 0;
   for (Long64_t i = 0; i < t->GetEntries() && done < maxTrk; ++i) {
      t->GetEntry(i);
      if (pe->GetEntries() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      for (auto &tr : ev->GetTrackCand()) {
         auto *hc = tr.GetHitClusterArray();
         const int n = hc->size();
         if (n < 20) continue;
         double s = 0, lo = 1e9, hi = -1e9;
         for (int k = 0; k < n; ++k) { double z = (*hc)[k].GetPosition().Z();
            lo = std::min(lo, z); hi = std::max(hi, z);
            if (k) s += std::fabs(z - (*hc)[k-1].GetPosition().Z()); }
         if (hi - lo < 20) continue;
         tv.push_back(s/(hi-lo)); ncl.push_back(n); ++done;
         if (done >= maxTrk) break;
      }
   }
   auto q=[&](double fr){ std::vector<double> v=tv; std::sort(v.begin(),v.end());
        return v.empty()?-1.0:v[std::min(v.size()-1,(size_t)(fr*v.size()))]; };
   long bad=0; for(double x:tv) if(x>1.5) ++bad;
   printf("\n  SIM 14C(d,p) 2.85 T, %zu tracks >= 20 clusters, z span > 20 mm\n", tv.size());
   printf("    z tv / z span:  p25 %.2f   median %.2f   p75 %.2f   p90 %.2f   max %.2f\n",
          q(0.25), q(0.50), q(0.75), q(0.90), q(1.0));
   printf("    frac > 1.5: %.1f %%\n", tv.empty()?0.0:100.0*bad/tv.size());
   printf("\n  by CLUSTER COUNT:\n    %-14s %6s %10s %12s\n","clusters","n","median tv","frac >1.5");
   const int clo[4]={20,100,300,600}, chi[4]={100,300,600,100000};
   for(int b=0;b<4;++b){ std::vector<double> s; long nb=0;
      for(size_t k=0;k<tv.size();++k) if(ncl[k]>=clo[b]&&ncl[k]<chi[b]){ s.push_back(tv[k]); if(tv[k]>1.5) ++nb; }
      if(s.size()<5) continue; std::sort(s.begin(),s.end());
      printf("    %5d - %-6d %6zu %10.2f %11.0f%%\n", clo[b], chi[b]>99999?9999:chi[b], s.size(),
             s[s.size()/2], 100.0*nb/s.size()); }
}
