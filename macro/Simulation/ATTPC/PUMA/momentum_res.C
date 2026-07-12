/// @file momentum_res.C
/// @brief Transverse-momentum resolution from the PRA circle fit (no Kalman fit):
///        p_T [GeV/c] = 0.3 * B[T] * R[m], B = 4 T. Compares reconstructed p_T to
///        the truth for the PUMA test-8 (pi) / test-10 (K) samples, baseline vs DLC.
/// Run: root -b -q 'momentum_res.C("data/reco_pi_base.root",0.3749,"pi base")'
#include <vector>
#include <algorithm>
double iqrSigma(std::vector<double> v){ if(v.size()<4) return 0; std::sort(v.begin(),v.end());
  return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double median(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end());
  return v[v.size()/2]; }
void momentum_res(TString file="data/reco_pi_base.root", double pTruth=0.3749, TString tag="pi base")
{
   const double B = 4.0; // Tesla
   TFile f(file); if(f.IsZombie()){printf("no file %s\n",file.Data());return;}
   auto*t=(TTree*)f.Get("cbmsim");
   TClonesArray*pat=nullptr; t->SetBranchAddress("AtPatternEvent",&pat);
   std::vector<double> pt, dpFrac;
   long long n=t->GetEntries();
   for(long long i=0;i<n;i++){ t->GetEntry(i);
      auto*pe=(AtPatternEvent*)(pat?pat->At(0):nullptr); if(!pe) continue;
      for(auto&trk: pe->GetTrackCand()){
         double R_mm=trk.GetGeoRadius(); if(R_mm<=0||R_mm>1e4) continue;
         double p=0.3*B*(R_mm/1000.0); // GeV/c
         if(p<0.02||p>2.0) continue;   // sane band
         pt.push_back(p); dpFrac.push_back((p-pTruth)/pTruth);
      }
   }
   printf(">>> PRES tag=%-8s tracks=%zu  median_pT=%.3f GeV  truth=%.3f  bias=%+.1f%%  sigma_pT/pT=%.1f%%\n",
      tag.Data(), pt.size(), median(pt), pTruth, 100*median(dpFrac), 100*iqrSigma(dpFrac));
}
