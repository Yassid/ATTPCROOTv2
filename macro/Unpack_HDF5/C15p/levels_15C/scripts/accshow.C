#include "binacc.h"
#include "binbuild.C"
#include "psmodel.h"
// What does the exact per-bin acceptance look like, and does dens*acc match the measured
// continuum shape inside each theta_CM bin?
void accshow(){
  TFile f("dd_cmbins.root");
  for(auto&k:CFG){
    Pars A=parsOf(k);
    Win w{14,86,k.keLo,k.keHi};
    printf("\n================ %s  (window theta_lab %.0f-%.0f, KE %.1f-%.1f) ================\n",
           k.nm,w.thLo,w.thHi,w.keLo,w.keHi);
    for(int i=0;i<5;i++){
      double lo=TCM0+i*TCMSTEP, hi=lo+TCMSTEP;
      TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",k.nm,i)); if(!h||h->GetEntries()<300) continue;
      printf(" bin %d  theta_CM %2.0f-%2.0f   N %6.0f\n",i,lo,hi,h->GetEntries());
      printf("    Ex :"); for(double X=0;X<=20.01;X+=2) printf("%7.0f",X); printf("\n");
      printf("    acc:"); for(double X=0;X<=20.01;X+=2) printf("%7.3f",accBin(A,k.eb,X,lo,hi,w)); printf("\n");
      printf("    dat:"); for(double X=0;X<=20.01;X+=2){ double s=0;int n=0;
        for(int b=1;b<=h->GetNbinsX();b++){double x=h->GetBinCenter(b);
          if(x>X-1&&x<X+1){s+=h->GetBinContent(b);n++;}} printf("%7.0f",n?s/n:0.0);} printf("\n");
    }
  }
}
