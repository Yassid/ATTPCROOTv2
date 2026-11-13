#include "binfit.C"
// How erratic are the per-bin amplitudes really?  Print the PEAK HEIGHTS (not significances)
// per theta_CM bin, plus their bin-to-bin ratios.  A real angular distribution varies smoothly;
// independent per-bin fits have nothing enforcing that.
void angdiag(){
  for(auto&k:CFG){
    USE_LEBKG=true; FREE_TURNON=true; TURNON=2.5;
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,9.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; USE_LEBKG=true; }
    printf("\n== %s   BW peak height per theta_CM bin (counts / 100 keV)\n",k.nm);
    printf("   bin  cm      N   ");
    for(int i=0;i<NB;i++) printf("  BW%d(%.2f)",i+1,EPP[i]-PPGS); printf("   chi2\n");
    double A[6][NB]={{0}};
    for(int ib=0;ib<5;ib++){
      BR o=fitBin(k,ib,-2.0,9.0,false,true);
      if(o.N<400) continue;
      printf("   %2d  %2.0f-%2.0f %6.0f   ",ib,o.cmLo,o.cmHi,o.N);
      for(int i=0;i<NB;i++){ A[ib][i]=o.amp[i]; printf(" %9.1f",o.amp[i]); }
      printf("  %5.2f\n",o.chi2ndf); }
    // smoothness: how much does each amplitude jump bin to bin, relative to its own mean?
    printf("   bin-to-bin |ratio| (a smooth angular distribution gives O(1) steps)\n   ");
    for(int i=0;i<NB;i++){
      double worst=0;
      for(int ib=0;ib<4;ib++){ double a=A[ib][i],b=A[ib+1][i];
        if(a>1&&b>1){ double r=a>b?a/b:b/a; worst=std::max(worst,r); }
        else if(a>5||b>5) worst=std::max(worst,99.0); }   // one of them is zero
      printf("  BW%d %5.1f",i+1,worst); }
    printf("\n");
    USE_LEBKG=false; FREE_TURNON=false;
  }
}
