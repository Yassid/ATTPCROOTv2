/// @file cmp_gf_dirs.C
/// @brief Compare two genfit productions track by track, matched on (entry, trackID).
///
/// Used as the NULL TEST for AtGenfitter::SetTruncatePercent: with truncPct = 0 the output must
/// be IDENTICAL to production. If it is not, no truncation result from this build means anything.
/// Also used to quantify what each truncation setting actually moves.
#include <cmath>
#include <map>
#include <vector>
void cmp_gf_dirs(TString run = "run_0020_multifit",
                 TString dirA = "/mnt/f/a1975/gf_dp_cateloss/",
                 TString dirB = "/mnt/f/a1975/gf_dp_trunc/null/",
                 TString lblA = "production", TString lblB = "truncPct=0")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto open=[&](TString d){ return TFile::Open(d + run + "_genfitter_p.root"); };
   auto *fa=open(dirA), *fb=open(dirB);
   if(!fa||fa->IsZombie()||!fb||fb->IsZombie()){ printf("  cannot open both\n"); return; }
   auto *ta=(TTree*)fa->Get("cbmsim"); auto *tb=(TTree*)fb->Get("cbmsim");
   TClonesArray *aa=nullptr,*ab=nullptr;
   ta->SetBranchAddress("AtTrackingEvent",&aa); tb->SetBranchAddress("AtTrackingEvent",&ab);
   Long64_t N=std::min(ta->GetEntries(), tb->GetEntries());
   printf("\n  %s: %lld entries   %s: %lld entries   comparing %lld\n",
          lblA.Data(), ta->GetEntries(), lblB.Data(), tb->GetEntries(), N);
   long matched=0, identical=0, onlyA=0, onlyB=0;
   std::vector<double> dke,dth,dc2, keA,keB;
   for(Long64_t i=0;i<N;++i){
      ta->GetEntry(i); tb->GetEntry(i);
      auto *ea=(AtTrackingEvent*)aa->At(0); auto *eb=(AtTrackingEvent*)ab->At(0);
      if(!ea||!eb) continue;
      std::map<int,AtFittedTrack*> ma, mb;
      for(auto &t: ea->GetFittedTracks()) ma[t->GetTrackID()]=t.get();
      for(auto &t: eb->GetFittedTracks()) mb[t->GetTrackID()]=t.get();
      for(auto &kv: ma){ if(!mb.count(kv.first)){ ++onlyA; continue; } }
      for(auto &kv: mb){ if(!ma.count(kv.first)){ ++onlyB; continue; } }
      for(auto &kv: ma){
         auto it=mb.find(kv.first); if(it==mb.end()) continue;
         ++matched;
         auto &ka=kv.second->GetKinematicsXtr(); auto &kb=it->second->GetKinematicsXtr();
         double a=ka.kineticEnergy, b=kb.kineticEnergy;
         double thA=ka.theta*TMath::RadToDeg(), thB=kb.theta*TMath::RadToDeg();
         double na=kv.second->GetTrackMetadata()->GetNdf(), nb=it->second->GetTrackMetadata()->GetNdf();
         double ca=na>0?kv.second->GetTrackMetadata()->GetChi2()/na:-1;
         double cb=nb>0?it->second->GetTrackMetadata()->GetChi2()/nb:-1;
         if(a==b && thA==thB) ++identical;
         dke.push_back(a-b); dth.push_back(thA-thB); dc2.push_back(ca-cb);
         keA.push_back(a); keB.push_back(b);
      }
   }
   auto med=[](std::vector<double> v){ if(v.empty())return 0.0; std::sort(v.begin(),v.end()); return v[v.size()/2]; };
   auto medabs=[&](std::vector<double> v){ for(auto &x:v) x=std::fabs(x); return med(v); };
   printf("  matched %ld   BIT-IDENTICAL (KE and theta) %ld = %.1f%%   only in %s %ld   only in %s %ld\n",
          matched, identical, matched?100.0*identical/matched:0.0, lblA.Data(), onlyA, lblB.Data(), onlyB);
   if(matched){
      printf("  median KE   %s %.4f   %s %.4f   median |dKE| %.5f MeV\n",
             lblA.Data(), med(keA), lblB.Data(), med(keB), medabs(dke));
      printf("  median |dtheta| %.5f deg    median d(chi2/ndf) %.4f\n", medabs(dth), med(dc2));
   }
   fa->Close(); fb->Close();
}
