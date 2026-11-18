#include "one1015.C"
// theta_CM 10-15, BOTH codes, TWO background cases, fitted g.s. to 10 MeV with no gap.
//   FREE        p3* q (Ex-Sn)^a exp(-k(Ex-Sn))   2 shape par -- the background can make a bump
//               near 4 MeV and absorb the 4.43/4.93 doublet
//   CONSTRAINED p3* q exp(-k(Ex-Sn))             1 shape par -- the turn-on is fixed by phase
//               space, so strength at 4-5 has to come from the STATES
// The question is whether P1 survives when the background cannot imitate it.
static double BP[2][2][16]; static int BN[2][2];
void two1015(int BIN=0);
void two1015(int BIN){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas(Form("c2_%d",BIN),"c2",1400,1000); c->Divide(2,2);
  const char* cn[2]={"FREE  (E_{x}-S_{n})^{a} e^{-k(E_{x}-S_{n})}","CONSTRAINED  e^{-k(E_{x}-S_{n})}"};
  int forms[2]={3,2}, nsh[2]={2,1};
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    B1::BIN=BIN; B1::EB=CFG[ic].eb;
    RK::Bin& d=RK::B[BIN];
    Pars A=parsOf(CFG[ic]); A.vzLo=300; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
    B1::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05;
      B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
    printf("\n================ %s  theta_CM %.0f-%.0f ================\n",CFG[ic].nm,d.cmLo,d.cmHi);
    printf("   %-14s %9s %9s %16s %9s %9s %9s %9s\n",
           "case","-2lnL","chi2/ndf","shape","Ags","P1","P2","7.4+8.3");
    double csol[16]; bool haveC=false;
    for(int ord=0;ord<2;ord++){
      int cs = ord==0 ? 1 : 0;              // constrained first, then free seeded from it
      B1::FORM=forms[cs]; B1::NSH=nsh[cs];
      int npar=3+nsh[cs]+4+(B1::GAUSS?2:0); BN[ic][cs]=npar; B1::NPAR=npar;
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(500000); m->SetTolerance(0.01); m->SetPrintLevel(0); m->SetStrategy(2);
      ROOT::Math::Functor fn(&B1::fcn,npar); m->SetFunction(fn);
      double pk=d.h->GetBinContent(d.h->GetMaximumBin());
      m->SetLimitedVariable(0,"Ags",pk,pk*0.01,0,pk*10);
      m->SetLimitedVariable(1,"A52",0.1*pk,pk*0.001,0,pk*5);
      m->SetLimitedVariable(2,"B",50,1,0,1e7);
      if(false){
        // free, seeded from the constrained solution: a starts at 0 (where the forms coincide)
        // ⚠ NOT 0.0: Minuit2's internal transformation is singular AT a limit, and seeding
        // exactly on the lower bound threw the ATTPCROOT fit to -2lnL 126260 (chi2 634).
        m->SetLimitedVariable(3,"s0",0.05,0.02,0.0,9.0);
        m->SetLimitedVariable(4,"s1",csol[3],0.02,0.0,9.0);
        for(int g=0;g<4;g++) m->SetLimitedVariable(5+g,Form("A%d",g),csol[4+g],
                                                   std::max(1e-3,csol[4+g]*0.01),0,1e5);
        m->SetVariableValue(2,csol[2]);
      } else {
        for(int s=0;s<nsh[cs];s++) m->SetLimitedVariable(3+s,Form("s%d",s),0.5,0.02,0.0,9.0);
        for(int g=0;g<4;g++) m->SetLimitedVariable(3+nsh[cs]+g,Form("A%d",g),20,0.5,0,1e5); }
      if(B1::GAUSS){
        m->SetLimitedVariable(npar-2,"mu",d.gsMu,0.005,d.gsMu-0.6,d.gsMu+0.6);
        m->SetLimitedVariable(npar-1,"sig",d.sgR,0.005,0.08,1.5); }
      m->Minimize();
      std::vector<double> bestp(m->X(),m->X()+npar); double bestL=m->MinValue();
      if(cs==0){
        // scan the extra shape parameter; a = 0 reproduces the constrained form exactly
        const double AG[8]={0.0,0.25,0.5,0.75,1.0,1.5,2.0,3.0};
        for(double av:AG){
          auto* ms=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
          ms->SetMaxFunctionCalls(500000); ms->SetTolerance(0.01); ms->SetPrintLevel(0);
          ms->SetStrategy(2);
          ROOT::Math::Functor fs(&B1::fcn,npar); ms->SetFunction(fs);
          double pk2=d.h->GetBinContent(d.h->GetMaximumBin());
          ms->SetLimitedVariable(0,"Ags",haveC?csol[0]:pk2,pk2*0.01,0,pk2*10);
          ms->SetLimitedVariable(1,"A52",haveC?csol[1]:0.1*pk2,pk2*0.001,0,pk2*5);
          ms->SetLimitedVariable(2,"B",haveC?csol[2]:50,1,0,1e7);
          ms->SetFixedVariable(3,"s0",av);
          ms->SetLimitedVariable(4,"s1",haveC?csol[3]:0.5,0.02,0.0,9.0);
          for(int g=0;g<4;g++) ms->SetLimitedVariable(5+g,Form("A%d",g),
              haveC?csol[4+g]:20.0,std::max(1e-3,(haveC?csol[4+g]:20.0)*0.01),0,1e5);
          if(B1::GAUSS){
            ms->SetLimitedVariable(npar-2,"mu",haveC?csol[npar-2]:d.gsMu,0.005,d.gsMu-0.6,d.gsMu+0.6);
            ms->SetLimitedVariable(npar-1,"sig",haveC?csol[npar-1]:d.sgR,0.005,0.08,1.5); }
          ms->Minimize();
          if(ms->MinValue()<bestL){ bestL=ms->MinValue();
            bestp.assign(ms->X(),ms->X()+npar); }
          delete ms; } }
      const double* p=bestp.data();
      for(int j=0;j<npar;j++) BP[ic][cs][j]=p[j];
      if(cs==1){ for(int j=0;j<npar;j++) csol[j]=p[j]; haveC=true; }
      double s2=0; int nd=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
        if(x<B1::exLo||x>B1::exHi) continue;
        double mm=B1::model(x,p); if(mm<1e-9) mm=1e-9; s2+=(n-mm)*(n-mm)/std::max(mm,1.0); nd++; }
      // peak heights of each group
      const double* a=p+3+nsh[cs];
      auto grpH=[&](int g){ double mx=0;
        for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==g&&d.in[r]){
          double mu=(RK::EPP6[r]-RK::PPGS6)+d.gsMu, sm=0;
          for(int q2=0;q2<RK::NR;q2++) if(RK::GRP[q2]==g&&d.in[q2])
            sm+=((RK::RIDX[q2]>=0)?RK::FIXRAT[RK::RIDX[q2]]:1.0)*d.T[q2].Eval(mu);
          mx=std::max(mx,std::max(0.0,a[g])*sm); }
        return mx; };
      double h0=grpH(0), h1=grpH(1), h23=grpH(2)+grpH(3);   // BEFORE delete m
      TString sh; for(int s=0;s<nsh[cs];s++) sh+=Form("%6.3f ",p[3+s]);
      printf("   %-14s %9.1f %9.2f %16s %9.0f %9.1f %9.1f %9.1f   mu %+.3f  FWHM %.3f\n",
             cs?"constrained":"free",bestL,s2/std::max(1.0,nd-(double)npar),
             sh.Data(),p[0],h0,h1,h23,B1::GAUSS?p[npar-2]:d.gsMu,
             2.355*(B1::GAUSS?p[npar-1]:d.sgR));
      delete m;
      // ---- draw
      c->cd(1+cs*2+ic); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone(Form("h%d%d%d",BIN,ic,cs)); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(B1::exLo,B1::exHi);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=B1::exLo&&x<=B1::exHi) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(4*std::max(1.0,mx)); hd->SetMinimum(0.4);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.55); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ   %s;E_{x} [MeV];counts / 100 keV",
                        CFG[ic].nm,d.cmLo,d.cmHi,cn[cs]));
      hd->Draw("E");
      int cols[7]={ic?kRed+1:kAzure+2,kViolet+1,kCyan+2,kGray+2,kGreen+3,kOrange+7,kMagenta+2};
      int stys[7]={1,1,1,7,1,1,1};
      for(int kk=0;kk<7;kk++){
        auto* g=new TGraph();
        for(double x=B1::exLo;x<=B1::exHi;x+=0.02){
          double Xt=x-(B1::GAUSS?BP[ic][cs][BN[ic][cs]-2]:d.gsMu), v=0; const double* aa=BP[ic][cs]+3+nsh[cs];
          if(kk==0) v=B1::model(x,BP[ic][cs]);
          else if(kk==1){ if(B1::GAUSS){
                            double mu=BP[ic][cs][BN[ic][cs]-2],sg=std::max(0.05,BP[ic][cs][BN[ic][cs]-1]);
                            double t=(x-mu)/sg; v=std::max(0.0,BP[ic][cs][0])*std::exp(-0.5*t*t); }
                          else v=std::max(0.0,BP[ic][cs][0])*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR); }
          else if(kk==2){ if(B1::GAUSS){
                            double mu=BP[ic][cs][BN[ic][cs]-2],sg=std::max(0.05,BP[ic][cs][BN[ic][cs]-1]);
                            double t=(x-mu-0.740)/sg; v=std::max(0.0,BP[ic][cs][1])*std::exp(-0.5*t*t); }
                          else v=std::max(0.0,BP[ic][cs][1])*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR); }
          else if(kk==3) v=std::max(0.0,BP[ic][cs][2])*B1::shp(Xt,BP[ic][cs]+3)
                           *(B1::ACC?std::max(0.0,B1::ACC->Eval(Xt)):1.0);
          else { int grp=kk-4;
                 for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==grp&&d.in[r])
                   v+=std::max(0.0,aa[grp])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(x);
                 if(grp==2) for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==3&&d.in[r])
                   v+=std::max(0.0,aa[3])*d.T[r].Eval(x); }
          if(v>0) g->SetPoint(g->GetN(),x,v); }
        if(g->GetN()<2) continue;
        g->SetLineColor(cols[kk]); g->SetLineStyle(stys[kk]); g->SetLineWidth(kk==0?3:2); g->Draw("L"); }
      hd->Draw("E SAME");
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.034);
      tx->DrawLatex(0.46,0.86,Form("#chi^{2}/ndf = %.2f",s2/std::max(1.0,nd-(double)npar)));
      tx->SetTextSize(0.030); tx->SetTextColor(kGreen+3);
      tx->DrawLatex(0.46,0.80,Form("P1 peak = %.1f",h0));
      tx->SetTextColor(kOrange+7);
      tx->DrawLatex(0.46,0.75,Form("P2 peak = %.1f",h1));
      { auto* l=new TLine(PS::SN+d.gsMu,0.4,PS::SN+d.gsMu,hd->GetMaximum()*0.3);
        l->SetLineStyle(3); l->SetLineColor(kGray+2); l->Draw(); }
    }
    CFG[0]=b0; CFG[1]=b1;
  }
  c->SaveAs(Form("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/bin%02.0f%02.0f_two_cases%s.png",TCM0+BIN*TCMSTEP,TCM0+(BIN+1)*TCMSTEP,B1::GAUSS?"_gauss":"_gt"));
  printf("\n   wrote bin%02.0f%02.0f_two_cases%s.png\n",TCM0+BIN*TCMSTEP,TCM0+(BIN+1)*TCMSTEP,B1::GAUSS?"_gauss":"_gt");
  VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
