/// @file ke_vs_geo_dp.C
/// @brief Score a genfit production against the GEOMETRY, track by track.
///
/// The pattern circle gives an energy that owes nothing to the fit: at field B a radius R means
/// p = 0.29979*B*R, so KE_geo = sqrt(p^2+m^2)-m. That is how run_0016 e11488 was caught -- circle
/// 117.9 mm -> 5.39 MeV against a fitted 166.8 MeV, a factor 31.
///
/// IMPORTANT -- this is a GROSS FAILURE DETECTOR, NOT A PRECISION METRIC. A single circle fitted
/// to a decelerating spiral averages a radius that shrinks along the track (0.964 -> 0.734 vertex
/// to stopping end, measured), so KE_geo is itself biased low by tens of percent for long tracks.
/// A ratio of 1.3 means nothing. A ratio of 31 means the fit diverged. So the quantity reported
/// is the CATASTROPHIC FRACTION -- |log10(KE_fit/KE_geo)| beyond a threshold -- which no plausible
/// circle bias can produce.
#include <cmath>
#include <vector>
void ke_vs_geo_dp(TString run="run_0016_multifit",
                  TString dirA="/mnt/f/a1975/gf_dp_cateloss/", TString lblA="baseline (z-sort)",
                  TString dirB="/mnt/f/a1975/gf_dp_clorder/",  TString lblB="cluster order",
                  double B=2.85, double mass=938.272, int minCl=20)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto scan=[&](TString dir, const char *lbl){
      auto *f=TFile::Open(dir+run+"_genfitter_p.root");
      if(!f||f->IsZombie()){ printf("  %-22s cannot open\n", lbl); return; }
      auto *t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
      t->SetBranchAddress("AtTrackingEvent",&a);
      std::vector<double> ratio; long nfit=0, noGeo=0;
      std::vector<std::pair<int,double>> byN;   // (nclusters, ratio)
      for(Long64_t i=0;i<t->GetEntries();++i){
         t->GetEntry(i); auto *ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
         std::map<int,double> Rgeo;   // trackID -> pattern radius
         std::map<int,int> Ncl;
         for(auto &tr: ev->GetTrackArray())
            if(tr.GetHitClusterArray()->size()>=(size_t)minCl) {
               Rgeo[tr.GetTrackID()]=tr.GetGeoRadius();
               Ncl[tr.GetTrackID()]=tr.GetHitClusterArray()->size(); }
         for(auto &ft: ev->GetFittedTracks()){
            auto it=Rgeo.find(ft->GetTrackID()); if(it==Rgeo.end()) continue;
            double R=it->second;
            if(!(R>5.0) || !std::isfinite(R)){ ++noGeo; continue; }
            double p=0.29979*B*R;                       // MeV/c with R in mm
            double keGeo=std::sqrt(p*p+mass*mass)-mass;
            double keFit=ft->GetKinematicsXtr().kineticEnergy;
            if(!(keFit>0) || !(keGeo>0)) { ++noGeo; continue; }
            ++nfit; ratio.push_back(keFit/keGeo);
            byN.emplace_back(Ncl[ft->GetTrackID()], keFit/keGeo);
         }
      }
      f->Close();
      if(ratio.empty()){ printf("  %-22s no tracks\n", lbl); return; }
      std::vector<double> s=ratio; std::sort(s.begin(),s.end());
      auto q=[&](double p){ return s[(size_t)(p*(s.size()-1))]; };
      long c2=0,c5=0,c10=0;
      for(double r:ratio){ if(r>2||r<0.5) ++c2; if(r>5||r<0.2) ++c5; if(r>10||r<0.1) ++c10; }
      printf("  %-22s n %6ld | median %6.2f  p10 %6.2f  p90 %7.2f | >2x or <0.5x %5.1f%% | >5x %5.1f%% | >10x %5.1f%%\n",
             lbl, nfit, q(0.5), q(0.1), q(0.9),
             100.0*c2/nfit, 100.0*c5/nfit, 100.0*c10/nfit);
      // SHORT tracks make few turns, so their single circle is not averaging a shrinking radius --
      // that is where KE_geo can be trusted and the ratio should sit near 1.
      const int lo[4]={20,50,100,300}, hi[4]={50,100,300,100000};
      for(int b=0;b<4;++b){
         std::vector<double> v;
         for(auto &pr:byN) if(pr.first>=lo[b] && pr.first<hi[b]) v.push_back(pr.second);
         if(v.size()<20) continue;
         std::sort(v.begin(),v.end());
         long f5=0; for(double r:v) if(r>5||r<0.2) ++f5;
         printf("      ncl %4d-%-6d n %6zu   median %6.2f   >5x %5.1f%%\n",
                lo[b], hi[b]>99999?9999:hi[b], v.size(), v[v.size()/2], 100.0*f5/v.size());
      }
   };
   printf("\n  %s  --  KE_fit / KE_geo (circle radius at B = %.2f T)\n\n", run.Data(), B);
   scan(dirA, lblA.Data());
   scan(dirB, lblB.Data());
   printf("\n  the last three columns are the CATASTROPHIC fractions -- no circle bias makes 5x\n");
}
