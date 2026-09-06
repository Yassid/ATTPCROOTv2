/// Is the z sort actually scrambling THIS track, or only making harmless local swaps?
/// Compares the stored cluster order with the z-sorted one: how far each cluster moves, the total
/// path each order traverses, and the per-cluster z advance against the drift quantum.
#include <cmath>
#include <vector>
#include <numeric>
void zorder_check_dp(Long64_t entry=11488,int tid=0,TString run="run_0016_multifit",
                     TString gfDir="/mnt/f/a1975/gf_dp_find/", double tbMM=1.84)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto*f=TFile::Open(gfDir+run+"_genfitter_p.root");
   auto*t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
   t->SetBranchAddress("AtTrackingEvent",&a); t->GetEntry(entry);
   auto*ev=(AtTrackingEvent*)a->At(0);
   for(auto &tr: ev->GetTrackArray()){
      if(tr.GetTrackID()!=tid) continue;
      auto*hc=tr.GetHitClusterArray(); const int n=hc->size();
      std::vector<ROOT::Math::XYZPoint> P; for(int i=0;i<n;++i) P.push_back((*hc)[i].GetPosition());
      std::vector<int> zi(n); std::iota(zi.begin(),zi.end(),0);
      std::sort(zi.begin(),zi.end(),[&](int x,int y){return P[x].Z()<P[y].Z();});
      // where does each cluster end up?
      std::vector<int> posInZ(n);
      for(int r=0;r<n;++r) posInZ[zi[r]]=r;
      long big=0; double maxMove=0;
      for(int i=0;i<n;++i){ double d=std::fabs((double)posInZ[i]-i);
         if(d>0.05*n) ++big; maxMove=std::max(maxMove,d); }
      auto pathOf=[&](const std::vector<int>&ord){ double s=0;
         for(size_t k=1;k<ord.size();++k) s+=(P[ord[k]]-P[ord[k-1]]).R(); return s; };
      std::vector<int> id(n); std::iota(id.begin(),id.end(),0);
      // per-cluster z advance along the STORED order
      std::vector<double> dz;
      for(int i=1;i<n;++i) dz.push_back(std::fabs(P[i].Z()-P[i-1].Z()));
      std::sort(dz.begin(),dz.end());
      // consecutive-in-z pairs: how far apart in 3D?
      std::vector<double> sep;
      for(int r=1;r<n;++r) sep.push_back((P[zi[r]]-P[zi[r-1]]).R());
      std::vector<double> ss=sep; std::sort(ss.begin(),ss.end());
      printf("\n  %s e%lld t%d : %d clusters, z span %.1f mm\n",run.Data(),entry,tid,n,
             P[zi[n-1]].Z()-P[zi[0]].Z());
      printf("\n  --- is z a usable sort key here? ---\n");
      printf("    mean z advance per cluster (stored order) : %.2f mm\n",
             (P[zi[n-1]].Z()-P[zi[0]].Z())/(n-1));
      printf("    median |dz| between consecutive stored    : %.2f mm   (one time bucket = %.2f)\n",
             dz[dz.size()/2], tbMM);
      printf("    fraction of steps with |dz| < 1 bucket    : %.0f%%\n",
             100.0*std::count_if(dz.begin(),dz.end(),[&](double d){return d<tbMM;})/dz.size());
      printf("\n  --- what the z sort actually does to the order ---\n");
      printf("    clusters moved by >5%% of the array        : %ld  (%.1f%%)\n", big, 100.0*big/n);
      printf("    largest move                              : %.0f positions (%.0f%% of the array)\n",
             maxMove, 100.0*maxMove/n);
      printf("    total path, stored order                  : %.0f mm\n", pathOf(id));
      printf("    total path, z-sorted order                : %.0f mm   (%.2fx)\n",
             pathOf(zi), pathOf(zi)/pathOf(id));
      printf("    3D gap between consecutive-in-z clusters  : median %.1f mm, p90 %.1f, max %.1f\n",
             ss[ss.size()/2], ss[(size_t)(0.9*ss.size())], ss.back());
   }
   f->Close();
}
