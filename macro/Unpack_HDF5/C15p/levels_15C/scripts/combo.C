#include "binfit.C"
// Combine the two fixes: restrict the upper limit (so the fit avoids the high-Ex tail with no
// physics) AND put the background turn-on ABOVE the g.s. region, where the turn-on scan wanted
// it anyway.  With the turn-on at Sn the background's rising edge sits under the g.s. tail, so
// improving Ex>2 costs Ex<2 -- that is what made bin 0 look so bad.
void chi2split2(const Cfg&k,int ib,double hi,double&c1,int&n1,double&c2r,int&n2,double&tot){
  BR o=fitBin(k,ib,-2.0,hi,false,true);
  c1=c2r=0; n1=n2=0; tot=o.chi2ndf; if(o.N<400||!o.h) return;
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
void combo(){
  for(auto&k:CFG){
    printf("\n================ %s ================\n",k.nm);
    printf("   hi  turn-on  bin   chi2(Ex<2)  chi2(Ex>2)   chi2/ndf   nBW\n");
    for(double hi : {7.0,8.0,9.0}){
      for(double to : {1.218,2.0,2.5}){
        USE_LEBKG=true; FREE_TURNON=(to!=1.218); TURNON=to;
        REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
        { double sa[NB]={0},tot=0; USE_LEBKG=false;
          for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,hi,false,false);
            if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
          for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; USE_LEBKG=true; }
        for(int ib=0;ib<2;ib++){
          double a,b2,tt; int na,nb; chi2split2(k,ib,hi,a,na,b2,nb,tt);
          if(!na&&!nb) continue;
          printf("  %4.1f   %5.2f    %d   %9.2f  %10.2f  %9.2f    %d\n",
                 hi,to,ib,na?a/na:0,nb?b2/nb:0,tt,NBW_IN); }
      }
      printf("\n");
    }
    USE_LEBKG=false; FREE_TURNON=false;
  }
}
