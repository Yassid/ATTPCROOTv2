#include "binfit.C"
// Does removing the low-vz background improve the per-bin fits, and by how much does it cut the
// Ex 2.5-4 excess?  Scan the low-vz cut.
void vzcutscan(){
  for(int ic=0;ic<2;ic++){
    Cfg&k=CFG[ic];
    printf("\n================ %s ================\n",k.nm);
    printf("  vzLo    N(bin0)   chi2/ndf bin0  bin1  bin2    excess frac (Ex2.5-4 / g.s.)\n");
    for(double v : {0.0,100.0,150.0,200.0,250.0,300.0}){
      VZLO=v;
      BB_QUIET=true; binbuild();   // rebuild the binned spectra with this cut
      REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
      { double sa[NB]={0},tot=0;
        for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
          if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
        for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
      double c[3]={0,0,0}; double n0=0,fr=0;
      for(int ib=0;ib<3;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,true);
        c[ib]=o.chi2ndf; if(ib==0){ n0=o.N;
          double g=o.h->Integral(o.h->FindBin(-0.5),o.h->FindBin(0.5));
          double e=o.h->Integral(o.h->FindBin(2.5),o.h->FindBin(4.0));
          fr=g>0?100.0*e/g:0; } }
      printf("  %5.0f  %8.0f   %8.2f %8.2f %8.2f      %6.2f %%\n",v,n0,c[0],c[1],c[2],fr);
    }
    VZLO=0;
  }
}
