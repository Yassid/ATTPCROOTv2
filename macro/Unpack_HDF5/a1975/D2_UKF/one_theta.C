#include "pra_theta.C"
void one_theta(Long64_t entry, int tid, TString run="run_0020_multifit",
               TString gfDir="/mnt/f/a1975/gf_dp_cateloss/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f=TFile::Open(gfDir+run+"_genfitter_p.root");
   auto *t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
   t->SetBranchAddress("AtTrackingEvent",&a); t->GetEntry(entry);
   auto *ev=(AtTrackingEvent*)a->At(0);
   for(auto &tr: ev->GetTrackArray()){
      if(tr.GetTrackID()!=tid) continue;
      auto o=praThetaHits(tr);
      printf("\n  entry %lld tid %d   hits %zu   clusters %zu\n",entry,tid,
             tr.GetHitArray().size(), tr.GetHitClusterArray()->size());
      printf("  stored GeoTheta   %.2f deg\n", tr.GetGeoTheta()*TMath::RadToDeg());
      if(!o.ok){printf("  replay: NO RANSAC LINE\n");return;}
      printf("  replay theta      %.2f deg   (slope sign %+d)\n", o.angleDeg,o.sign);
      printf("  circle            R %.1f mm  centre (%.1f, %.1f)   consensus %d hits\n",
             o.radius,o.cx,o.cy,o.nUsedCircle);
      printf("  arc-vs-z points   %zu   RANSAC inliers %d  (%.0f%%)\n", o.arc.size(),
             (int)std::count(o.inlier.begin(),o.inlier.end(),1),
             100.0*std::count(o.inlier.begin(),o.inlier.end(),1)/o.arc.size());
      double alo=*std::min_element(o.arc.begin(),o.arc.end()), ahi=*std::max_element(o.arc.begin(),o.arc.end());
      printf("  arc range         [%.0f, %.0f] mm   line dir (%.3f, %.3f)\n",alo,ahi,o.dirX,o.dirY);
   }
}
