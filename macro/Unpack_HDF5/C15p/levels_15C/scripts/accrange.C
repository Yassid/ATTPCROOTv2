#include "seqps.C"
// acceptance(Ex) with each code's OWN (theta, KE) window, not the 40-90 placeholder.
void accrange(){
  struct C{const char* nm; double thLo,thHi,keLo,keHi,kb;};
  C c[2]={{"ATTPCROOT",14,86,0.5,18,191.25},{"Spyral",14,86,0.0,20,195.0}};
  for(auto&k:c){
    ACC_TH_LO=k.thLo; ACC_TH_HI=k.thHi; ACC_KE_LO=k.keLo; ACC_KE_HI=k.keHi;
    printf("\n%-10s theta %.0f-%.0f  KE %.1f-%.1f\n",k.nm,k.thLo,k.thHi,k.keLo,k.keHi);
    printf("   Ex :"); for(double X=2;X<=22.01;X+=2) printf("%7.0f",X); printf("\n");
    printf("   acc:"); for(double X=2;X<=22.01;X+=2) printf("%7.3f",acc_of(X,k.kb)); printf("\n");
    printf("   d*a:");
    double n9=PS::dens(9,k.kb)*acc_of(9,k.kb);
    for(double X=2;X<=22.01;X+=2) printf("%7.2f",PS::dens(X,k.kb)*acc_of(X,k.kb)/n9); printf("   (norm at 9)\n");
  }
  // the background level is MEASURABLE where (d,d') cannot reach: high Ex
  TFile f("dd_ex.root");
  for(const char* nm:{"ATTPCROOT","Spyral"}){ TH1D* h=(TH1D*)f.Get(nm);
    auto avg=[&](double a,double b){ double s=0;int n=0;
      for(int i=1;i<=h->GetNbinsX();i++){double x=h->GetBinCenter(i);
        if(x>a&&x<b){s+=h->GetBinContent(i);n++;}} return n?s/n:0.0; };
    printf("\n%-10s data: Ex 14-18 %.1f/bin   18-22 %.1f/bin   below Sn %.1f/bin\n",
           nm,avg(14,18),avg(18,22),avg(-2.5,-1.5)); }
}
