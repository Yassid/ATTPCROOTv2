#include "ridge.C"
// Verify the 4th band's identity properly. t, 3He and 4He rails have DIFFERENT SHAPES in Brho,
// so measuring the band's position at several Brho and comparing residuals against each candidate
// discriminates -- a wrong species diverges even if it matches at one point.
void verify(){
  TH2F* h2=loadPlane("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root",1045,1225,1.0,"vf");
  int ib=h2->GetYaxis()->FindBin(0.60);
  TH1D* px=h2->ProjectionX("vd",ib,ib); px->SetDirectory(nullptr);
  int bm=0; double vm=-1;
  for(int b=px->FindBin(2.0);b<=px->FindBin(12.0);++b) if(px->GetBinContent(b)>vm){vm=px->GetBinContent(b);bm=b;}
  double K0=px->GetBinCenter(bm)/sqrt(betheRel(0.60,1,1));
  printf("K0 = %.4f from proton %.2f at Brho 0.60\n\n", K0, px->GetBinCenter(bm));
  struct S{const char*n;double A,z;};
  S cand[4]={{"t",3,1},{"3He",3,2},{"4He",4,2},{"6He",6,2}};
  // follow the 4th band empirically: seed near the 4He prediction at Brho 0.60, then free-run with
  // a tight window, and record where it actually goes.
  double seed=K0*sqrt(betheRel(0.60,4,2));
  printf("seeding the 4th band at %.2f (the 4He prediction) and following it freely\n",seed);
  printf("%8s %9s | %s\n","Brho","observed","residual (obs-pred)/pred for each candidate");
  printf("%8s %9s |","","");
  for(auto&c:cand) printf(" %9s",c.n); printf("\n");
  double g=seed; int nrow=0; double sres[4]={0,0,0,0}; int ns=0;
  for(int i=h2->GetYaxis()->FindBin(0.60); i<=h2->GetNbinsY(); ++i){
    double br=h2->GetYaxis()->GetBinCenter(i); if(br>1.25) break;
    Pt q=fitBin(h2,i,g,0.10); if(!q.ok) continue;
    g=q.mu;
    if(nrow%4==0){
      printf("%8.3f %9.2f |",br,q.mu);
      for(int k=0;k<4;k++){ double pr=K0*sqrt(betheRel(br,cand[k].A,cand[k].z));
        printf(" %+8.1f%%",100*(q.mu-pr)/pr); }
      printf("\n");
    }
    for(int k=0;k<4;k++){ double pr=K0*sqrt(betheRel(br,cand[k].A,cand[k].z));
      sres[k]+=fabs((q.mu-pr)/pr); }
    ns++; nrow++;
  }
  printf("\nmean |residual| over %d bins:",ns);
  for(int k=0;k<4;k++) printf("  %s %.1f%%",cand[k].n,100*sres[k]/std::max(1,ns));
  printf("\n-> the band belongs to the species with the smallest mean residual\n");
}
