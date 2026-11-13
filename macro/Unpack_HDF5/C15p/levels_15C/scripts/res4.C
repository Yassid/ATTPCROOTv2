#include "resfit.C"
// Final configuration: Ex 3-9, 5 deg bins to 60, the five populated resonances with the
// unresolvable 4.43/4.93 doublet TIED.  Four angular distributions: D12, 5.98, 6.68, 7.38.
void res4(){
  setCM(10,5,60);
  RF::USE[0]=RF::USE[1]=RF::USE[2]=RF::USE[3]=RF::USE[4]=true; RF::USE[5]=false;
  RF::TIE12=true;
  for(auto&k:CFG){
    std::vector<std::vector<double>> amp;
    rfRun(k,10,true,&amp);
    const char* nm[4]={"D12 4.4+4.9","5.98       ","6.68       ","7.38       "};
    printf("   solid-angle corrected, each normalised to its own max:\n");
    printf("   theta_CM   :");
    for(int b=0;b<(int)amp.size();b++) printf(" %6.0f",TCM0+(b+0.5)*TCMSTEP); printf("\n");
    for(int j=0;j<4;j++){
      std::vector<double> v; double mx=0;
      for(int b=0;b<(int)amp.size();b++){
        double s=std::sin((TCM0+(b+0.5)*TCMSTEP)*M_PI/180.0);
        double x=amp[b][j]/s; v.push_back(x); mx=std::max(mx,x); }
      printf("   %s:",nm[j]);
      for(double x:v) printf(" %6.3f",mx>0?x/mx:0); printf("\n"); }
  }
  RF::TIE12=false;
}
