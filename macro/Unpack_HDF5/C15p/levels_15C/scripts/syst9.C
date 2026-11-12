#include "binfit.C"
// The background systematic AGAIN, in the adopted window [-2, 9].  With the turn-on moved to 2.5
// the background now sits directly under the resonance region, so this has to be re-checked:
// a better chi2 is worth nothing if it is bought by making the resonances arbitrary.
void syst9(){
  for(auto&k:CFG){
    printf("\n================ %s  (theta_CM bin 0 and 1, Ex [-2,9]) ================\n",k.nm);
    for(int ib=0;ib<2;ib++){
      REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
      { double sa[NB]={0},tot=0; USE_LEBKG=false; USE_LONGTAIL=false;
        for(int j=0;j<6;j++){ BR o=fitBin(k,j,-2.0,9.0,false,false);
          if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
        for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
      printf("   bin %d   model                chi2  |  BW peak heights\n",ib);
      double mn[NB],mx[NB]; for(int i=0;i<NB;i++){mn[i]=1e30;mx[i]=-1e30;}
      auto rec=[&](const char* nm,BR&o,bool keep){
        printf("           %-20s %6.2f |",nm,o.chi2ndf);
        for(int i=0;i<NB;i++){ printf(" %7.1f",o.amp[i]);
          if(keep){ mn[i]=std::min(mn[i],o.amp[i]); mx[i]=std::max(mx[i],o.amp[i]); } }
        printf("%s\n",keep?"":"   (not counted)"); };
      USE_LEBKG=false; USE_LONGTAIL=false; FREE_TURNON=false;
      { BR o=fitBin(k,ib,-2.0,9.0,false,true); rec("none",o,true); }
      USE_LEBKG=true;
      for(double to : {1.218,2.0,2.5,3.0}){ FREE_TURNON=(to!=1.218); TURNON=to;
        BR o=fitBin(k,ib,-2.0,9.0,false,true);
        bool keep=(o.bgK<1.0&&o.chi2ndf<6);
        rec(Form("erf @ %.2f",to),o,keep); }
      FREE_TURNON=false; USE_LEBKG=false; USE_LONGTAIL=true;
      { BR o=fitBin(k,ib,-2.0,9.0,false,true); rec("long elastic tail",o,true); }
      USE_LONGTAIL=false;
      printf("           %-20s        | spread/mean:","");
      for(int i=0;i<NB;i++){ double m=0.5*(mn[i]+mx[i]);
        printf(" %6.0f%%",(m>1&&mx[i]>-1e29)?100.0*(mx[i]-mn[i])/m:0.0); }
      printf("\n\n"); }
  }
}
