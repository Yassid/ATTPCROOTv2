/// @file fitcmp_order_dp.C
/// @brief Judge cluster orderings by THE FIT, not by geometry proxies.
///
/// Geometry metrics are biased here: total path length is minimised BY CONSTRUCTION by a greedy
/// nearest-neighbour walk, so it cannot fairly compare against AtPRA's walk. The honest referee is
/// the outcome -- feed each ordering to the SAME fitter and score the fitted energy against the
/// one the pattern circle implies (p = 0.29979*B*R), which owes nothing to the fit.
///
/// Reports the GROSS FAILURE fraction (>5x off) per ordering. On long tracks that threshold is not
/// reachable by circle-fit noise; on short ones it is, so this is run on long tracks only.
#include "order_helix_dp.C"
#include "vertex_order_fit_dp.C"
#include <fstream>
void fitcmp_order_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                     TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=120,
                     int minCl=100, double B=2.85, double mass=938.272,
                     TString geoName="ATTPC_D300torr_v2_geomanager.root",
                     TString eLossTable="proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf=TFile::Open(dir+"/geometry/"+geoName); gf->Get("FAIRGeom"); }
   TString eloss = dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable;
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   std::vector<double> rS, rZ, rH;   // KE_fit / KE_geo per ordering
   long okS=0, okZ=0, okH=0, refused=0;
   for(auto &j:jobs){
      if(done>=maxTracks) break;
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if(rt!=open){ if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); a=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&a); open=rt; }
      t->GetEntry(en); auto*ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
      for(auto &tr: ev->GetTrackArray()){
         if(tr.GetTrackID()!=id) continue;
         auto*hc=tr.GetHitClusterArray(); const int n=hc->size(); if(n<minCl) continue;
         const std::vector<AtHitCluster> stored=*hc;
         std::vector<ROOT::Math::XYZPoint> P; P.reserve(n);
         for(int i=0;i<n;++i) P.push_back(stored[i].GetPosition());
         double R=0,h=0,c=0; auto ord=orderHelix(P,h,&R,&c);
         if(!(R>5)) continue;
         double p=0.29979*B*R, keGeo=std::sqrt(p*p+mass*mass)-mass;
         if(!(keGeo>0)) continue;
         ++done;
         auto score=[&](const std::vector<AtHitCluster>&u, std::vector<double>&dst, long &ok){
            VoOut o=voFit(*&tr,u,eloss);
            if(o.ok && o.ke>0){ dst.push_back(o.ke/keGeo); ++ok; } };
         score(stored, rS, okS);
         { std::vector<AtHitCluster> u=stored;
           std::sort(u.begin(),u.end(),[](const AtHitCluster&A,const AtHitCluster&B){
              return A.GetPosition().Z()<B.GetPosition().Z(); });
           score(u, rZ, okZ); }
         if(ord.empty()) ++refused;
         else { std::vector<AtHitCluster> u; u.reserve(ord.size());
                for(int i:ord) u.push_back(stored[i]); score(u, rH, okH); }
         if(done%20==0) printf("  %d tracks...\n", done);
      }
   }
   if(f) f->Close();
   auto rep=[&](const char*lbl, std::vector<double>&v, long ok){
      if(v.empty()){ printf("  %-20s no fits\n", lbl); return; }
      std::vector<double> s=v; std::sort(s.begin(),s.end());
      long f5=0,f2=0; for(double x:v){ if(x>5||x<0.2) ++f5; if(x>2||x<0.5) ++f2; }
      printf("  %-20s fits %4ld   median %6.2f   >2x %5.1f%%   >5x %5.1f%%\n",
             lbl, ok, s[s.size()/2], 100.0*f2/v.size(), 100.0*f5/v.size()); };
   printf("\n  %d long tracks (>=%d clusters), scored against the circle energy\n", done, minCl);
   printf("  helix refused on %ld\n\n", refused);
   rep("stored (PRA walk)", rS, okS);
   rep("z sort", rZ, okZ);
   rep("HELIX PHASE", rH, okH);
}
