/// Compare three cluster orderings on the same tracks: the stored PRA walk, a z sort, and the
/// intrinsic helix phase. Judged on three order-free quantities:
///   * HOPS  -- steps beyond 3x the mean step. A correct order has none.
///   * TOTAL PATH -- a shorter total path means a shorter traversal.
///   * Z TOTAL VARIATION / Z SPAN -- 1.0 for a monotonic traversal, 2.0 for one full U-turn.
///
/// !!! TOTAL PATH IS A BIASED REFEREE HERE AND MUST NOT DECIDE ANYTHING ON ITS OWN.  The stored
/// PRA order IS a greedy nearest-neighbour walk, i.e. a path-length minimiser, so scoring it on
/// path length scores it on its own objective and it wins by construction.  Worse, the failure
/// this file exists to find -- the walk cutting across turns of a spiral and U-turning in z -- is
/// a way of achieving a SHORT path with a WRONG order, so path length cannot see it at all.  The
/// header of the first version asserted "the true trajectory is the shortest sensible traversal",
/// and on a multi-turn spiral that is simply false.
///
/// The z total variation is the referee that is not the greedy walk's objective: z is monotonic
/// along any real helix, so any excess over 1.0 is either z noise or a wrong order, and comparing
/// orderings on the SAME tracks holds the noise fixed.
#include "order_helix_dp.C"
#include <fstream>
/// z monotonicity: total variation of z along the order, over the z span.  Order-sensitive and
/// NOT what the greedy walk optimises, which is the whole point of having it.
static double zTV(const std::vector<ROOT::Math::XYZPoint> &p)
{
   const int n = p.size(); if (n < 3) return -1;
   double tv = 0, lo = 1e9, hi = -1e9;
   for (int i = 0; i < n; ++i) { double z = p[i].Z(); lo = std::min(lo,z); hi = std::max(hi,z);
      if (i) tv += std::fabs(z - p[i-1].Z()); }
   if (hi - lo < 20) return -1;           // no z lever arm; the ratio is noise over noise
   return tv / (hi - lo);
}
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
   // PAIRED: only tracks the helix ordering accepted go into both, so the two columns
   // describe the same tracks and the comparison is not a change of sample.
   std::vector<double> tvS, tvH;
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
         double vS=zTV(P);
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
         // the z referee, on the SAME tracks so the z noise is held fixed.  A z SORT is
         // monotonic by construction and is therefore not scored: 1.000 would mean nothing.
         double vH=zTV(byH);
         if(vS>0 && vH>0){ tvS.push_back(vS); tvH.push_back(vH); }
      }
   }
   if(f) f->Close();
   auto q=[](std::vector<double> v,double f){ if(v.empty()) return -1.0; std::sort(v.begin(),v.end());
        return v[std::min(v.size()-1,(size_t)(f*v.size()))]; };
   auto frac=[](const std::vector<double>&v,double thr){ if(v.empty()) return -1.0; long n=0;
        for(double x:v) if(x>thr) ++n; return 100.0*n/v.size(); };
   auto med=[](std::vector<double> v){ if(v.empty()) return -1.0; std::sort(v.begin(),v.end()); return v[v.size()/2]; };
   printf("\n  %d tracks  (helix REFUSED on %d = %.1f%%: no pitch reached concentration 0.5)\n\n",
          done, refused, done?100.0*refused/done:0.0);
   printf("  %-22s %10s %12s %14s\n","ordering","tracks w/ hop","total hops","median path / stored");
   printf("  %-22s %10ld %12ld %14s\n","stored (PRA walk)",trS,hopS,"1.000");
   printf("  %-22s %10ld %12ld %14.3f\n","z sort",trZ,hopZ,med(ratZ));
   printf("  %-22s %10ld %12ld %14.3f\n","HELIX PHASE",trH,hopH,med(ratH));
   printf("\n  median pitch concentration %.3f  (1 = the points lie exactly on one helix)\n", med(conc));
   printf("  PATH LENGTH IS A BIASED REFEREE -- the stored order is itself a path minimiser, and\n"
          "  cutting across turns is a way to get a SHORT path with a WRONG order.  Read the z\n"
          "  total variation instead:\n\n");
   printf("  %-22s %8s %8s %8s %8s %10s\n","z tv / z span","n","p25","median","p75","frac >1.5");
   printf("  %-22s %8zu %8.2f %8.2f %8.2f %9.0f%%\n","stored (PRA walk)",tvS.size(),
          q(tvS,0.25),q(tvS,0.50),q(tvS,0.75),frac(tvS,1.5));
   printf("  %-22s %8zu %8.2f %8.2f %8.2f %9.0f%%\n","HELIX PHASE",tvH.size(),
          q(tvH,0.25),q(tvH,0.50),q(tvH,0.75),frac(tvH,1.5));
   printf("\n  1.0 is a monotonic traversal, 2.0 one full U-turn.  The helix column is also the\n"
          "  NOISE FLOOR of the metric: it is the same clusters with the same z errors, so\n"
          "  whatever it sits at is what a correct order costs, and the excess over it in the\n"
          "  stored column is ordering error and nothing else.\n");
}
