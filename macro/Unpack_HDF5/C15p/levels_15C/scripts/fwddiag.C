#include "resfitC.C"
// Why does the 6.68 yield collapse in the forward two theta_CM bins, which hold the MOST data?
// Suspect: the ground-state tail.  Those bins have the largest elastic peak, and the response is
// a bifurcated Gaussian with EXPONENTIAL tails that reach all the way under Ex 3-9.  If the tail
// is overestimated it eats the resonance.  Compare, at Ex = 6.68, the three model pieces against
// the data.
void fwddiag(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALON; RH::TIE12=true;
  bool keep[6]={true,true,true,true,false,false};
  for(int i=0;i<RH::NR;i++) RH::USE[i]=keep[i];
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
    std::vector<std::vector<double>> out;
    rhRun(CFG[ic],NCM,false,&out);
    printf("\n======== %s : anatomy of the model at Ex = 6.68 ========\n",CFG[ic].nm);
    printf("   %-9s %8s %9s  | %8s %8s %8s %8s | %8s %8s\n",
           "theta_CM","N","data@6.7","contin.","D12 tail","5.98 tail","6.68 BW","model","resid");
    for(int b=0;b<NCM;b++){
      RH::Bin& d=RH::B[b]; if(!d.h) continue;
      double X=(6.68-0.119)+d.gsMu;  // where the 6.68 peak sits
      int ib=d.h->FindBin(X); double dat=d.h->GetBinContent(ib);
      double Xt=X-d.gsMu;
      double gs = d.gsAmp*shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
      double e5 = d.e5Amp*shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
      double bg = RH::LASTBG[b]*RH::bshape(Xt,RH::LASTP0)*(d.ga?d.ga->Eval(Xt):1.0);
      double a443=d.in[0]?RH::LASTA[b][0]*d.T[0].Eval(X):0.0;   // the TIED doublet, at 6.68
      double a493=d.in[1]?RH::LASTA[b][0]*d.T[1].Eval(X):0.0;
      double a598=d.in[2]?RH::LASTA[b][2]*d.T[2].Eval(X):0.0;   // 5.98 BW tail reaching 6.68
      double bw  =d.in[3]?RH::LASTA[b][3]*d.T[3].Eval(X):0.0;   // the 6.68 itself
      double tot=gs+e5+bg+a443+a493+a598+bw;
      printf("   %3.0f-%-5.0f %8.0f %9.1f  | %8.2f %8.2f %8.2f %8.2f | %8.2f %8.2f\n",
             d.cmLo,d.cmHi,d.N,dat,bg,a443+a493,a598,bw,tot,dat-tot); }
    printf("   response of the forward bins: ");
    for(int b=0;b<2;b++) printf("bin%d FWHM_L %.2f FWHM_R %.2f alL %.2f alR %.2f gsAmp %.0f   ",
        b,2.355*RH::B[b].sgL,2.355*RH::B[b].sgR,RH::B[b].alL,RH::B[b].alR,RH::B[b].gsAmp);
    printf("\n");
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true; BB_QUIET=true; binbuild();
}
