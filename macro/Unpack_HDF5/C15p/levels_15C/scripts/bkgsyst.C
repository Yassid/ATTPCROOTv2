#include "binfit.C"
// How much do the resonance amplitudes move when the low-energy background model is varied?
// Scan the turn-on centre over the range where the fit is STABLE (background still a background:
// B bounded, kappa < 0.4, resonances not collapsed) and record the spread in each BW.
void bkgsyst(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    printf("\n================ %s  (theta_CM bin 0) ================\n",k.nm);
    printf("   background model          chi2   |  BW peak heights (4.43 4.93 5.98 6.68 7.38 8.28)\n");
    double mn[NB],mx[NB]; for(int i=0;i<NB;i++){mn[i]=1e30;mx[i]=-1e30;}
    auto rec=[&](const char* nm,BR&o){
      printf("   %-24s %6.2f  |",nm,o.chi2ndf);
      for(int i=0;i<NB;i++){ printf(" %7.1f",o.amp[i]);
        mn[i]=std::min(mn[i],o.amp[i]); mx[i]=std::max(mx[i],o.amp[i]); }
      printf("\n"); };
    USE_LEBKG=false; USE_LONGTAIL=false; { BR o=fitBin(k,0,-2.0,12.0,false,true); rec("none",o); }
    USE_LEBKG=true;  FREE_TURNON=false;  { BR o=fitBin(k,0,-2.0,12.0,false,true); rec("erf at Sn = 1.218",o); }
    FREE_TURNON=true;
    for(double c : {0.0,1.0,2.0,2.5}){ TURNON=c;
      BR o=fitBin(k,0,-2.0,12.0,false,true);
      if(o.bgK>0.45||o.chi2ndf>6) { printf("   %-24s  (unstable, skipped)\n",Form("erf at %.1f",c)); continue; }
      rec(Form("erf at %.1f",c),o); }
    FREE_TURNON=false; USE_LEBKG=false; USE_LONGTAIL=true;
    { BR o=fitBin(k,0,-2.0,12.0,false,true); rec("long elastic tail",o); }
    USE_LONGTAIL=false;
    printf("   ------------------------------------------------------------------\n");
    printf("   spread / mean            %%      |");
    for(int i=0;i<NB;i++){ double m=0.5*(mn[i]+mx[i]);
      printf(" %6.0f%%",m>0?100.0*(mx[i]-mn[i])/m:0.0); }
    printf("\n");
  }
}
