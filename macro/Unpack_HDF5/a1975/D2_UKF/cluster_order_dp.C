/// @file cluster_order_dp.C
/// @brief Is the stored CLUSTER array actually contiguous along the track?
///
/// AtPRA::OrderClustersAlongTrack (AtPRA.cxx:338) is a GREEDY nearest-neighbour walk seeded at the
/// highest-Z cluster. On a tightening spiral, successive turns come closer together than the
/// along-track step, so the walk can HOP TO THE NEXT TURN and come back later. When that happens,
/// "the first N clusters" is not a contiguous arc from the vertex -- which is why truncating can
/// leave a piece of the far end alive.
///
/// Prints per cluster: the 3D step from the previous one, and the unwrapped azimuth/arc, so a hop
/// shows up as a big step AND a jump in arc.
///   root -b -q 'cluster_order_dp.C(604,0)'
#include <cmath>
#include <vector>
void cluster_order_dp(Long64_t entry, int tid, TString run = "run_0020_multifit",
                      TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/", double jumpFac = 3.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f = TFile::Open(gfDir + run + "_genfitter_p.root");
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *a = nullptr; t->SetBranchAddress("AtTrackingEvent", &a);
   t->GetEntry(entry);
   auto *ev = (AtTrackingEvent *)a->At(0); if (!ev) { printf("no event\n"); return; }
   for (auto &tr : ev->GetTrackArray()) {
      if (tr.GetTrackID() != tid) continue;
      auto *hc = tr.GetHitClusterArray();
      const int n = hc->size();
      if (n < 5) { printf("only %d clusters\n", n); return; }
      // circle from the clusters, to get an arc coordinate independent of the array order
      double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
      for (int i=0;i<n;++i){auto p=(*hc)[i].GetPosition();double x=p.X(),y=p.Y();
         Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y;}
      double C=n*Sxx-Sx*Sx,D=n*Sxy-Sx*Sy,E=n*Sxxx+n*Sxyy-(Sxx+Syy)*Sx;
      double G=n*Syy-Sy*Sy,H=n*Sxxy+n*Syyy-(Sxx+Syy)*Sy, den=2*(C*G-D*D);
      double cx=0,cy=0,R=0;
      if (std::fabs(den)>1e-9){ cx=(E*G-D*H)/den; cy=(C*H-D*E)/den;
         for(int i=0;i<n;++i){auto p=(*hc)[i].GetPosition(); R+=std::hypot(p.X()-cx,p.Y()-cy);} R/=n; }
      std::vector<double> step(n,0), phi(n,0);
      double meanStep=0;
      for (int i=0;i<n;++i){
         auto p=(*hc)[i].GetPosition();
         phi[i]=std::atan2(p.Y()-cy,p.X()-cx)*TMath::RadToDeg();
         if(i){ auto q=(*hc)[i-1].GetPosition(); step[i]=(p-q).R(); meanStep+=step[i]; }
      }
      meanStep/= (n-1);
      printf("\n  %s entry %lld tid %d : %d clusters,  circle R %.1f mm centre (%.1f,%.1f)\n",
             run.Data(), entry, tid, n, R, cx, cy);
      printf("  mean step along the stored order: %.1f mm   (a HOP is flagged at > %.1fx that)\n\n",
             meanStep, jumpFac);
      printf("  %5s %8s %8s %8s %8s %9s %9s\n","idx","x","y","z","|xy|","step[mm]","phi[deg]");
      int nJump=0;
      for (int i=0;i<n;++i){
         auto p=(*hc)[i].GetPosition();
         bool hop = (i>0 && step[i] > jumpFac*meanStep);
         if (hop) ++nJump;
         if (i<8 || i>n-8 || hop)
            printf("  %5d %8.1f %8.1f %8.1f %8.1f %9.1f %9.1f %s\n", i, p.X(), p.Y(), p.Z(),
                   std::hypot(p.X(),p.Y()), step[i], phi[i], hop?"  <== HOP":"");
      }
      printf("\n  %d hops (> %.1fx mean step) in the stored cluster order\n", nJump, jumpFac);
      // does the first half of the array cover a contiguous arc?
      double zf=1e9,zl=-1e9, zf2=1e9,zl2=-1e9;
      for(int i=0;i<n/2;++i){double z=(*hc)[i].GetPosition().Z(); zf=std::min(zf,z); zl=std::max(zl,z);}
      for(int i=n/2;i<n;++i){double z=(*hc)[i].GetPosition().Z(); zf2=std::min(zf2,z); zl2=std::max(zl2,z);}
      printf("  z span of FIRST half of the array : [%.1f, %.1f]\n", zf, zl);
      printf("  z span of SECOND half             : [%.1f, %.1f]   %s\n", zf2, zl2,
             (zf2 < zl && zf < zl2) ? "<== THE TWO HALVES OVERLAP IN z" : "");
   }
}
