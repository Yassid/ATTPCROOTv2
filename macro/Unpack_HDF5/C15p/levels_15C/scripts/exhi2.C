#include "binfit.C"
// (a) where does each bin's chi2 come from -- the g.s. peak or the resonance region?
// (b) the upper-limit scan with the low-energy background switched ON.
void chi2split(const Cfg&k,int ib,double hi,double&c1,int&n1,double&c2r,int&n2){
  BR o=fitBin(k,ib,-2.0,hi,false,true);
  c1=c2r=0; n1=n2=0; if(o.N<400||!o.h) return;
  double fl=o.floor,cs=o.cs,gsA=o.gsA,e5A=o.e5A,gmu=GSMU,sl=SGL,sr=SGR,al=ALF_L,ar=ALF_R;
  double bB=o.bgB,bW=o.bgW,bK=o.bgK;
  std::vector<T2> T(NB); for(int i=0;i<NB;i++) T[i]=TT[i];
  TGraph* GCl=GC; bool in[NB]; for(int i=0;i<NB;i++) in[i]=BW_IN[i];
  auto mdl=[&](double X){ double v=fl+cs*std::max(0.0,GCl->Eval(X))
      +gsA*shapeGT(X,gmu,sl,sr,al,ar)+e5A*shapeGT(X,gmu+E52d,sl,sr,al,ar);
    if(USE_LEBKG){ double c0=FREE_TURNON?TURNON:B2::SN;
      v+=std::max(0.0,bB)*0.5*(1+std::erf((X-c0)/(M_SQRT2*std::max(0.05,bW))))*std::exp(-std::max(0.0,bK)*X); }
    for(int i=0;i<NB;i++) if(in[i]) v+=o.A[i]*T[i].Eval(X); return v; };
  for(int b=1;b<=o.h->GetNbinsX();b++){ double x=o.h->GetBinCenter(b),y=o.h->GetBinContent(b);
    if(x<-2.0||x>hi||y<10) continue;
    double d=(y-mdl(x))*(y-mdl(x))/y;
    if(x<2.0){c1+=d;n1++;} else {c2r+=d;n2++;} }
}
void exhi2(){
  for(auto&k:CFG){
    printf("\n================ %s ================\n",k.nm);
    for(int mode=0;mode<2;mode++){
      USE_LEBKG=(mode==1);
      printf("  low-energy background %s\n",mode?"ON ":"OFF");
      printf("   hi   bin    chi2/bin  Ex<2 (g.s.)   chi2/bin  Ex>2 (resonances)\n");
      for(double hi : {7.0,8.0,9.0,12.0}){
        REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
        { double sa[NB]={0},tot=0; bool sv=USE_LEBKG; USE_LEBKG=false;
          for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,hi,false,false);
            if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
          for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; USE_LEBKG=sv; }
        for(int ib=0;ib<3;ib++){
          double a,b2; int na,nb; chi2split(k,ib,hi,a,na,b2,nb);
          if(!na&&!nb) continue;
          printf("  %4.1f   %d     %8.2f (%2d bins)      %8.2f (%2d bins)\n",
                 hi,ib,na?a/na:0,na,nb?b2/nb:0,nb); }
      }
    }
    USE_LEBKG=false;
  }
}
