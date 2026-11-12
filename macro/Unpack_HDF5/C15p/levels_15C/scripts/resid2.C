#include "binfit.C"
// What is the SHAPE of what the model misses?  Fit as now (tailed response + 6 BW + sequential
// phase space + fixed floor), then look at data - model bin by bin in Ex.
void resid2(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    printf("\n================ %s : data - model ================\n",k.nm);
    printf("  bin   Ex:");
    for(double X=1.5;X<=11.01;X+=0.5) printf("%6.1f",X); printf("\n");
    for(int ib=0;ib<3;ib++){
      BR o=fitBin(k,ib,-2.0,12.0,false,true);
      if(o.N<400) continue;
      double fl=o.floor,cs=o.cs,gsA=o.gsA,e5A=o.e5A,gmu=GSMU,sl=SGL,sr=SGR,al=ALF_L,ar=ALF_R;
      std::vector<T2> T(NB); for(int i=0;i<NB;i++) T[i]=TT[i];
      TGraph* GCl=GC;
      auto mdl=[&](double X){ double v=fl+cs*std::max(0.0,GCl->Eval(X))
        +gsA*shapeGT(X,gmu,sl,sr,al,ar)+e5A*shapeGT(X,gmu+E52d,sl,sr,al,ar);
        for(int i=0;i<NB;i++) v+=o.A[i]*T[i].Eval(X); return v; };
      TH1D* h=o.h;
      printf("  %2d  dat :",ib);
      for(double X=1.5;X<=11.01;X+=0.5) printf("%6.0f",h->GetBinContent(h->FindBin(X))); printf("\n");
      printf("      mod :");
      for(double X=1.5;X<=11.01;X+=0.5) printf("%6.0f",mdl(X)); printf("\n");
      printf("      resid:");
      for(double X=1.5;X<=11.01;X+=0.5){ double d=h->GetBinContent(h->FindBin(X))-mdl(X);
        printf("%6.0f",d); } printf("\n");
      printf("      pull :");
      for(double X=1.5;X<=11.01;X+=0.5){ double y=h->GetBinContent(h->FindBin(X));
        double d=y-mdl(X); printf("%6.1f",y>0?d/std::sqrt(y):0.0); } printf("\n\n");
    }
  }
}
