#include "binfit.C"
// Re-check: with the LOW-ENERGY BACKGROUND now in the model, does the "missing strength at
// Ex ~ 3.4" survive?  If the background absorbs it, it was never a level.
static T2 TX3;
double fBinX3(double* x, double* p){ return fBin(x,p) + std::max(0.0,p[12])*TX3.Eval(x[0]); }
void extrabw2(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    printf("\n================ %s ================\n",k.nm);
    for(int mode=0;mode<2;mode++){
      USE_LEBKG=(mode==1);
      printf("  low-energy background %s\n",mode? "ON ":"OFF");
      double base[2];
      for(int ib=0;ib<2;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,true); base[ib]=o.nll; }
      printf("    Ex_extra:"); for(double E=2.0;E<=6.01;E+=0.4) printf("%8.1f",E); printf("\n");
      printf("    d(-2lnL):");
      for(double E=2.0;E<=6.01;E+=0.4){
        double sum=0;
        for(int ib=0;ib<2;ib++){
          BR o=fitBin(k,ib,-2.0,12.0,false,true);
          TX3.build(E+GSMU,0.20,1,SGR,-3.0,23.0);
          TFile f("dd_cmbins.root"); TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",k.nm,ib));
          h=(TH1D*)h->Clone(Form("q%s%d%d%d",k.nm,ib,(int)(E*10),mode)); h->SetDirectory(nullptr);
          TF1 g(Form("gq%s%d%d%d",k.nm,ib,(int)(E*10),mode),fBinX3,-2.0,12.0,13);
          g.SetParameter(0,o.gsA); g.SetParLimits(0,0,1e9);
          g.SetParameter(1,o.e5A); g.SetParLimits(1,0,1e9);
          g.SetParameter(2,o.cs);  g.SetParLimits(2,0,1e7);
          for(int i=0;i<NB;i++){ double s2=REF[i]*o.N/REFTOT;
            g.SetParameter(3+i,s2); g.SetParLimits(3+i,0,std::max(3.0*s2,0.4*REF[i])); }
          double tx=TX3.Eval(E+GSMU);
          // multi-start, as in fitBin: a single seed gives minima WORSE than the nested model
          // (d(-2lnL) of -3400 appeared here before this was added).
          double bestv=1e30;
          int nseed = USE_LEBKG ? 6 : 1;
          for(int sd=0; sd<nseed; sd++){
            g.SetParameter(0,o.gsA); g.SetParLimits(0,0,1e9);
            g.SetParameter(1,o.e5A); g.SetParLimits(1,0,1e9);
            g.SetParameter(2,o.cs);  g.SetParLimits(2,0,1e7);
            for(int i=0;i<NB;i++){ double s2=REF[i]*o.N/REFTOT;
              g.SetParameter(3+i,s2); g.SetParLimits(3+i,0,std::max(3.0*s2,0.4*REF[i])); }
            if(USE_LEBKG){
              const double FB[3]={0.3,1.0,2.0}, FK[2]={0.5,1.5};
              g.SetParameter(9, std::max(1.0,o.bgB)*FB[sd%3]); g.SetParLimits(9,0,1e6);
              g.SetParameter(10,std::max(0.2,o.bgW));          g.SetParLimits(10,0.05,3.0);
              g.SetParameter(11,std::max(0.02,o.bgK)*FK[sd/3]); g.SetParLimits(11,0.0,2.0);
            } else { g.FixParameter(9,0); g.FixParameter(10,1); g.FixParameter(11,0); }
            g.SetParameter(12,tx>0?10.0/tx:1); g.SetParLimits(12,0,1e9);
            TFitResultPtr r=h->Fit(&g,"RQ0NLS");
            double v=2*r->MinFcnValue(); if(v<bestv) bestv=v; }
          sum += base[ib]-bestv; }
        printf("%8.1f",sum); }
      printf("\n");
    }
    USE_LEBKG=false;
  }
}
