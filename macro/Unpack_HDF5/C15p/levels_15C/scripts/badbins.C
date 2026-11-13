#include "resfit.C"
// What is wrong in the two bins with bad chi2?
void badbins(){
  setCM(10,5,60);
  RF::USE[0]=RF::USE[1]=RF::USE[2]=RF::USE[3]=RF::USE[4]=true; RF::USE[5]=false;
  struct {int ic,b;} bad[2]={{0,8},{1,9}};
  for(auto&z:bad){
    Cfg& k=CFG[z.ic];
    RF::B.clear(); for(int ib=0;ib<10;ib++) rfPrep(k,ib);
    RF::Bin& d=RF::B[z.b];
    printf("\n== %s bin %d (theta_CM %.0f-%.0f), N = %.0f\n",k.nm,z.b,d.cmLo,d.cmHi,d.N);
    printf("   Ex :"); for(double X=2.5;X<=9.51;X+=0.5) printf("%6.1f",X); printf("\n");
    printf("   N  :"); for(double X=2.5;X<=9.51;X+=0.5) printf("%6.0f",d.h->GetBinContent(d.h->FindBin(X))); printf("\n");
    // total counts in the window, and how many bins have <5
    double tot=0; int nlow=0,n=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i);
      if(x<3||x>9) continue; n++; tot+=d.h->GetBinContent(i);
      if(d.h->GetBinContent(i)<5) nlow++; }
    printf("   counts in Ex 3-9: %.0f over %d bins (mean %.1f/bin); %d bins have < 5 counts\n",
           tot,n,n?tot/n:0,nlow);
    printf("   g.s. tail at Ex=3: %.2f   (gsAmp %.0f, FWHM_R %.3f, alpha_R %.2f)\n",
           d.gsAmp*shapeGT(3.0,d.gsMu,d.sgL,d.sgR,d.alL,d.alR),d.gsAmp,2.355*d.sgR,d.alR);
  }
}
