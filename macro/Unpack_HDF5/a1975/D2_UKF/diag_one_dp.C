/// Everything the reconstruction knows about ONE track, from a genfit production.
#include <cmath>
void diag_one_dp(Long64_t entry, int tid, TString run="run_0016_multifit",
                 TString gfDir="/mnt/f/a1975/gf_dp_find/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f=TFile::Open(gfDir+run+"_genfitter_p.root");
   auto *t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
   t->SetBranchAddress("AtTrackingEvent",&a); t->GetEntry(entry);
   auto *ev=(AtTrackingEvent*)a->At(0); if(!ev){printf("no event\n");return;}
   printf("\n=== %s entry %lld ===  tracks in event: %zu   fitted: %zu\n",
          run.Data(), entry, ev->GetTrackArray().size(), ev->GetFittedTracks().size());
   for (auto &tr : ev->GetTrackArray()){
      auto *hc=tr.GetHitClusterArray();
      double zlo=1e9,zhi=-1e9,xymin=1e9,xymax=-1e9,qsum=0;
      for (size_t i=0;i<hc->size();++i){ auto p=(*hc)[i].GetPosition();
         zlo=std::min(zlo,p.Z()); zhi=std::max(zhi,p.Z());
         double xy=std::hypot(p.X(),p.Y()); xymin=std::min(xymin,xy); xymax=std::max(xymax,xy);
         qsum+=(*hc)[i].GetCharge(); }
      printf("\n  track %d: %zu hits, %zu clusters   GeoTheta %7.2f deg  R %6.1f mm  charge %.3g\n",
             tr.GetTrackID(), tr.GetHitArray().size(), hc->size(),
             tr.GetGeoTheta()*TMath::RadToDeg(), tr.GetGeoRadius(), qsum);
      printf("           z [%.1f, %.1f]   |xy| [%.1f, %.1f]\n", zlo,zhi,xymin,xymax);
      // ordering hops
      if (hc->size()>2){
         double mean=0; for(size_t i=1;i<hc->size();++i) mean+=((*hc)[i].GetPosition()-(*hc)[i-1].GetPosition()).R();
         mean/=(hc->size()-1); int nh=0; double mx=0;
         for(size_t i=1;i<hc->size();++i){ double d=((*hc)[i].GetPosition()-(*hc)[i-1].GetPosition()).R();
            if(d>3*mean){++nh; mx=std::max(mx,d);} }
         printf("           ordering: mean step %.1f mm, %d hops, largest %.1f mm\n", mean,nh,mx);
         auto c0=(*hc)[0].GetPosition(), cN=(*hc)[hc->size()-1].GetPosition();
         printf("           cl[0]  z %.1f |xy| %.1f      cl[last] z %.1f |xy| %.1f\n",
                c0.Z(),std::hypot(c0.X(),c0.Y()), cN.Z(),std::hypot(cN.X(),cN.Y()));
      }
   }
   for (auto &ft : ev->GetFittedTracks()){
      auto &k=ft->GetKinematicsXtr(); auto &kf=ft->GetKinematics();
      double ndf=ft->GetTrackMetadata()->GetNdf(), c2=ft->GetTrackMetadata()->GetChi2();
      auto v=ft->GetVertex();
      printf("\n  FIT tid %d: KE %8.3f MeV (raw %8.3f)  theta %7.2f deg (raw %7.2f)\n",
             ft->GetTrackID(), k.kineticEnergy, kf.kineticEnergy,
             k.theta*TMath::RadToDeg(), kf.theta*TMath::RadToDeg());
      printf("           chi2 %.1f / ndf %.0f = %.3f   vertex (%.1f, %.1f, %.1f) |xy| %.1f mm\n",
             c2, ndf, ndf>0?c2/ndf:-1, v.X(),v.Y(),v.Z(), std::hypot(v.X(),v.Y()));
   }
   f->Close();
}
