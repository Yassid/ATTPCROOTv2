#include "ridge.C"
// IDENTIFY every band by asking which (A,z) rail it sits on. The dE/dx scale is calibrated once
// from the proton, so for any candidate species the locus is fully predicted with no free
// parameter. Count 3-sigma detections along each candidate rail: the real species light up.
void whichband(){
  struct S{const char*n;double A,z;};
  S cand[14]={{"p",1,1},{"d",2,1},{"t",3,1},{"3He",3,2},{"4He",4,2},{"6He",6,2},
              {"6Li",6,3},{"7Li",7,3},{"7Be",7,4},{"9Be",9,4},{"10B",10,5},{"11B",11,5},
              {"12C",12,6},{"15C",15,6}};
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  for(auto&p:ps){
    TH2F* h2=loadPlane(p.f,p.il,p.ih,p.sc,Form("wb_%s",p.n)); if(!h2) continue;
    double sb0=0.60; int ib=h2->GetYaxis()->FindBin(sb0);
    TH1D* px=h2->ProjectionX("wbd",ib,ib); px->SetDirectory(nullptr);
    int bm=0; double vm=-1;
    for(int b=px->FindBin(2.0);b<=px->FindBin(12.0);++b) if(px->GetBinContent(b)>vm){vm=px->GetBinContent(b);bm=b;}
    double muP=px->GetBinCenter(bm);
    double K0 = muP / sqrt(betheRel(sb0,1,1));     // universal dE/dx calibration from the proton
    printf("\n===== %s   proton anchor %.2f at Brho %.2f  ->  K0 = %.4f\n", p.n, muP, sb0, K0);
    printf("   %-6s %8s %7s %7s %9s %9s\n","species","mu@0.60","bins","3sig","<S/B>","<sig/mu>");
    for(auto&c:cand){
      double mu0 = K0*sqrt(betheRel(sb0,c.A,c.z));
      if(mu0<1.0||mu0>44.0){ printf("   %-6s %8.2f   (off scale)\n",c.n,mu0); continue; }
      auto tr=followRail(h2,sb0,mu0,c.A,0.11);
      // NB followRail's internal rail uses z=1; redo the rail properly for z!=1
      int n3=0, nall=0; double ssb=0, sw=0;
      for(int ibn=1; ibn<=h2->GetNbinsY(); ++ibn){
        double br=h2->GetYaxis()->GetBinCenter(ibn); if(br<0.30||br>1.10) continue;
        double pred=K0*sqrt(betheRel(br,c.A,c.z)); if(pred<1.0||pred>44.0) continue;
        Pt q=fitBin(h2,ibn,pred,0.11); if(!q.ok) continue;
        if(fabs(q.mu-pred)/pred > 0.09) continue;
        nall++;
        double s3=q.amp/std::max(1.0,sqrt(std::max(q.bkg,1.0)));
        if(s3>3.0 && q.sg/q.mu>0.015 && q.sg/q.mu<0.25){ n3++; sw+=q.sg/q.mu; if(q.bkg>0) ssb+=q.amp/q.bkg; }
      }
      printf("   %-6s %8.2f %7d %7d %9.2f %9.4f\n", c.n, mu0, nall, n3, n3?ssb/n3:0, n3?sw/n3:0);
    }
    delete h2;
  }
}
