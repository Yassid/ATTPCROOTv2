#include "ridge.C"
// Detection on a real criterion: the band counts as SEEN in a Brho bin only if the Gaussian
// amplitude exceeds 3x the Poisson fluctuation of the local background (amp > 3*sqrt(bkg)), and
// the fitted width is physical. "A fit converged" is not a detection.
void run_ridge3(){
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  const char* bn[3]={"proton","deuteron","triton"}; const double A[3]={1,2,3};
  printf("%-10s %-9s | %5s %5s | %-17s | %-15s\n","plane","band","bins","3sig","<sigma/mu> +- rms","<significance>");
  double W[2][3], Wr[2][3]; int NB[2][3];
  for(int ip=0;ip<2;ip++){ auto&p=ps[ip];
    TH2F* h2=loadPlane(p.f,p.il,p.ih,p.sc,Form("r3_%s",p.n)); if(!h2) continue;
    double sb0=0.60; int ib=h2->GetYaxis()->FindBin(sb0);
    TH1D* px=h2->ProjectionX("s3d",ib,ib); px->SetDirectory(nullptr);
    int bm=0; double vm=-1;
    for(int b=px->FindBin(2.0);b<=px->FindBin(12.0);++b) if(px->GetBinContent(b)>vm){vm=px->GetBinContent(b);bm=b;}
    double muP=px->GetBinCenter(bm), rp=betheRel(sb0,1,1);
    for(int s=0;s<3;s++){
      double mu0=muP*sqrt(betheRel(sb0,A[s],1)/rp);
      auto tr=followRail(h2,sb0,mu0,A[s],0.13);
      std::vector<double> w, sig;
      int nall=0;
      for(auto&q:tr){ if(q.brho<0.30||q.brho>1.10) continue; nall++;
        double s3 = q.amp/std::max(1.0,sqrt(std::max(q.bkg,1.0)));
        if(s3>3.0 && q.sg/q.mu>0.015 && q.sg/q.mu<0.25){ w.push_back(q.sg/q.mu); sig.push_back(s3); } }
      double mw=0,rw=0,ms=0;
      for(double v:w) mw+=v; if(w.size()) mw/=w.size();
      for(double v:w) rw+=(v-mw)*(v-mw); if(w.size()>1) rw=sqrt(rw/(w.size()-1));
      for(double v:sig) ms+=v; if(sig.size()) ms/=sig.size();
      W[ip][s]=mw; Wr[ip][s]=rw; NB[ip][s]=w.size();
      printf("%-10s %-9s | %5d %5zu | %7.4f +- %-6.4f | %15.1f\n",
             p.n,bn[s],nall,w.size(), mw, rw, ms);
    }
    delete h2;
  }
  printf("\n=== deuteron width, ATTPCROOT vs Spyral ===\n");
  double d=W[0][1]-W[1][1];
  double e=sqrt(Wr[0][1]*Wr[0][1]/std::max(1,NB[0][1]) + Wr[1][1]*Wr[1][1]/std::max(1,NB[1][1]));
  printf("  %.4f vs %.4f  ->  difference %+.4f +- %.4f  (%.1f sigma)\n", W[0][1], W[1][1], d, e, e>0?fabs(d)/e:0);
}
