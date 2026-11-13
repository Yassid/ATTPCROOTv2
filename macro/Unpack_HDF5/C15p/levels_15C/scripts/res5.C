#include "resfit.C"
// The five resonances the data populates (8.28 dropped: gain 15.5 for Spyral, below the ~18
// threshold, and 24.9 for ATTPCROOT -- it is the one resonance the two codes disagree about).
// Angular distributions on the 10-point, 5 deg grid.
void res5(){
  setCM(10,5,60);
  RF::USE[0]=RF::USE[1]=RF::USE[2]=RF::USE[3]=RF::USE[4]=true; RF::USE[5]=false;
  for(auto&k:CFG){
    std::vector<std::vector<double>> amp;
    rfRun(k,10,true,&amp);
    printf("   solid-angle corrected (amp / sin theta), each column normalised to its own max:\n");
    printf("   theta_CM :");
    for(int b=0;b<(int)amp.size();b++) printf(" %6.0f",TCM0+(b+0.5)*TCMSTEP); printf("\n");
    for(int j=0;j<5;j++){
      std::vector<double> v; double mx=0;
      for(int b=0;b<(int)amp.size();b++){
        double s=std::sin((TCM0+(b+0.5)*TCMSTEP)*M_PI/180.0);
        double x=amp[b][j]/s; v.push_back(x); mx=std::max(mx,x); }
      printf("   %5.2f MeV:",RF::EPP6[j]-RF::PPGS6);
      for(double x:v) printf(" %6.3f",mx>0?x/mx:0); printf("\n"); }
  }
}
