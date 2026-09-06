/// Compare three cluster orderings on the same tracks: the stored PRA walk, a z sort, and the
/// intrinsic helix phase. Judged on two order-free quantities:
///   * HOPS  -- steps beyond 3x the mean step. A correct order has none.
///   * TOTAL PATH -- the true trajectory is the shortest sensible traversal, so a shorter total
///     path means a better ordering. This is the metric the earlier work should have used: the
///     median step is blind to a minority of scrambled pairs.
#include "order_helix_dp.C"
#include <fstream>
static void pathStats(const std::vector<ROOT::Math::XYZPoint>&p, double &total, int &nHop, double &maxHop)
{
   total=0; nHop=0; maxHop=0;
   const int n=p.size(); if(n<2) return;
   std::vector<double> s; s.reserve(n-1);
   for(int i=1;i<n;++i){ double d=(p[i]-p[i-1]).R(); s.push_back(d); total+=d; }
   double mean=total/s.size();
   for(double d:s) if(d>3*mean){ ++nHop; maxHop=std::max(maxHop,d); }
}
void test_order_helix_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                         TString gfDir="/mnt/f/a1975/gf_dp_cateloss/", int maxTracks=400)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0, refused=0;
   long hopS=0,hopZ=0,hopH=0; long trS=0,trZ=0,trH=0;
   std::vector<double> ratZ, ratH, conc;
   for(auto &j:jobs){
      if(done>=maxTracks) break;
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if(rt!=open){ if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); a=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&a); open=rt; }
      t->GetEntry(en); auto*ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
      for(auto &tr: ev->GetTrackArray()){
         if(tr.GetTrackID()!=id) continue;
         auto*hc=tr.GetHitClusterArray(); const int n=hc->size(); if(n<20) continue;
         std::vector<ROOT::Math::XYZPoint> P; P.reserve(n);
         for(int i=0;i<n;++i) P.push_back((*hc)[i].GetPosition());
         ++done;
         double tS,tZ,tH,mS,mZ,mH; int hS,hZ,hH;
         pathStats(P,tS,hS,mS);
         auto byZ=P; std::sort(byZ.begin(),byZ.end(),
            [](const ROOT::Math::XYZPoint&A,const ROOT::Math::XYZPoint&B){return A.Z()<B.Z();});
         pathStats(byZ,tZ,hZ,mZ);
         double h=0,R=0,c=0; auto ord=orderHelix(P,h,&R,&c);
         conc.push_back(c);
         if(ord.empty()){ ++refused; if(hS) ++trS; hopS+=hS; ratZ.push_back(tZ/tS); continue; }
         std::vector<ROOT::Math::XYZPoint> byH; byH.reserve(n);
         for(int i:ord) byH.push_back(P[i]);
         pathStats(byH,tH,hH,mH);
         hopS+=hS; hopZ+=hZ; hopH+=hH;
         if(hS) ++trS; if(hZ) ++trZ; if(hH) ++trH;
         if(tS>0){ ratZ.push_back(tZ/tS); ratH.push_back(tH/tS); }
      }
   }
   if(f) f->Close();
   auto med=[](std::vector<double> v){ if(v.empty()) return -1.0; std::sort(v.begin(),v.end()); return v[v.size()/2]; };
   printf("\n  %d tracks  (helix REFUSED on %d = %.1f%%: no pitch reached concentration 0.5)\n\n",
          done, refused, done?100.0*refused/done:0.0);
   printf("  %-22s %10s %12s %14s\n","ordering","tracks w/ hop","total hops","median path / stored");
   printf("  %-22s %10ld %12ld %14s\n","stored (PRA walk)",trS,hopS,"1.000");
   printf("  %-22s %10ld %12ld %14.3f\n","z sort",trZ,hopZ,med(ratZ));
   printf("  %-22s %10ld %12ld %14.3f\n","HELIX PHASE",trH,hopH,med(ratH));
   printf("\n  median pitch concentration %.3f  (1 = the points lie exactly on one helix)\n", med(conc));
   printf("  a SHORTER path than stored means a better ordering; hops should go to ZERO.\n");
}
