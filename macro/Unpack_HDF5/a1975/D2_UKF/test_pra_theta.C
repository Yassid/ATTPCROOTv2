/// Headless check of pra_theta.C: does the replay reproduce the stored GeoTheta?
#include "pra_theta.C"
void test_pra_theta(TString run="run_0020_multifit", TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int nMax=12)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f = TFile::Open(gfDir + run + "_genfitter_p.root");
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *a = nullptr; t->SetBranchAddress("AtTrackingEvent", &a);
   printf("\n  %8s %4s %6s %10s %10s %8s %10s\n","entry","tid","nhits","stored","replay","sign","inliers");
   int n=0;
   for (Long64_t i=0; i<t->GetEntries() && n<nMax; ++i) {
      t->GetEntry(i);
      auto *ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
      for (auto &tr : ev->GetTrackArray()) {
         if (tr.GetHitArray().size() < 100) continue;
         auto o = praThetaHits(tr);
         printf("  %8lld %4d %6zu %10.2f %10.2f %8d %6d/%zu\n", i, tr.GetTrackID(),
                tr.GetHitArray().size(), tr.GetGeoTheta()*TMath::RadToDeg(),
                o.ok?o.angleDeg:-1.0, o.sign, (int)std::count(o.inlier.begin(),o.inlier.end(),1), o.arc.size());
         if (++n>=nMax) break;
      }
   }
}
