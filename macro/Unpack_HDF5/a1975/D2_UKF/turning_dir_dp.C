/// @file turning_dir_dp.C
/// @brief Certify the vertex end from the SENSE OF ROTATION, which is fixed by q and B.
///
/// Yassid's point: for a known charge in a known field the trajectory must turn one specific way.
/// So the sign of the cumulative azimuth swept ALONG THE STORED ARRAY says whether the array runs
/// in the direction of travel or backwards -- with no seed, no |xy| heuristic and no ordering
/// assumption. AtPRA already uses this relation, but BACKWARDS: AtPRA.cxx:184-199 infers the
/// CHARGE from the rotation sense, ordering hits by distance from the origin as a "time proxy" --
/// an assumption a curling spiral breaks. Here q is known (proton) and the direction is inferred.
///
/// SELF-CALIBRATING: the sense is a property of q*B, so it is IDENTICAL for every proton in the
/// run. Whatever sign the majority shows IS the physical sense; the minority are the tracks whose
/// stored array runs backwards. That sidesteps the sim/data field-sign convention entirely
/// (bFieldSign +1 sim, -1 data) rather than assuming it -- though the convention still has to be
/// pinned before this is used to SET anything.
#include <cmath>
#include <vector>
#include <fstream>
void turning_dir_dp(TString listFile="/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                    TString gfDir="/mnt/f/a1975/gf_dp_find/", int maxTracks=800)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   TString open; TFile*f=nullptr; TTree*t=nullptr; TClonesArray*a=nullptr; int done=0;
   // [sense<0 ?][theta>90 ?]
   long tab[2][2]={{0,0},{0,0}};
   // does the sense agree with "the end nearer the beam axis is the vertex"?
   long agree=0, disagree=0, ambig=0;
   std::vector<double> turns;
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
         double cx=tr.GetGeoCenter().first, cy=tr.GetGeoCenter().second;
         if(!std::isfinite(cx)||!std::isfinite(cy)) continue;
         // cumulative unwrapped azimuth along the STORED array
         double cum=0, prev=0; bool first=true;
         for(int i=0;i<n;++i){ auto p=(*hc)[i].GetPosition();
            double ph=std::atan2(p.Y()-cy,p.X()-cx);
            if(!first){ double d=ph-prev; while(d>M_PI)d-=2*M_PI; while(d<=-M_PI)d+=2*M_PI; cum+=d; }
            prev=ph; first=false; }
         if(std::fabs(cum) < 0.5) continue;      // too little rotation to have a sense
         ++done;
         turns.push_back(cum/(2*M_PI));
         const int neg = (cum<0)?1:0;
         const double gt = tr.GetGeoTheta()*TMath::RadToDeg();
         const int bwd = (std::isfinite(gt) && gt>90)?1:0;
         tab[neg][bwd]++;
         // cross-check against the axis marker, on tracks where that marker is decisive
         auto c0=(*hc)[0].GetPosition(), cN=(*hc)[n-1].GetPosition();
         double xy0=std::hypot(c0.X(),c0.Y()), xyN=std::hypot(cN.X(),cN.Y());
         double lo=std::min(xy0,xyN), hi=std::max(xy0,xyN);
         if(lo>1e-6 && hi/lo>2.0){ if((xy0<xyN)==(neg==1)) ++agree; else ++disagree; }
         else ++ambig;
      }
   }
   if(f) f->Close();
   long nneg=tab[1][0]+tab[1][1], npos=tab[0][0]+tab[0][1];
   printf("\n  %d tracks with a measurable rotation sense\n\n", done);
   printf("  cumulative azimuth along the STORED array:\n");
   printf("    NEGATIVE (clockwise seen from +z) : %6ld  (%.1f%%)   of which backward %ld\n",
          nneg, 100.0*nneg/done, tab[1][1]);
   printf("    POSITIVE                          : %6ld  (%.1f%%)   of which backward %ld\n",
          npos, 100.0*npos/done, tab[0][1]);
   printf("\n  => the MAJORITY sense is the physical one for a proton in this field;\n");
   printf("     the %ld minority tracks (%.1f%%) have their array stored BACKWARDS.\n",
          std::min(nneg,npos), 100.0*std::min(nneg,npos)/done);
   std::sort(turns.begin(),turns.end());
   printf("\n  turns swept: p10 %.2f  median %.2f  p90 %.2f\n",
          turns[turns.size()/10], turns[turns.size()/2], turns[9*turns.size()/10]);
   printf("\n  cross-check against 'the end nearer the beam axis is the vertex' (decisive tracks):\n");
   printf("    agree %ld   disagree %ld   -> %.0f%% consistent   (%ld ambiguous excluded)\n",
          agree, disagree, (agree+disagree)? 100.0*agree/(agree+disagree):0.0, ambig);
}
