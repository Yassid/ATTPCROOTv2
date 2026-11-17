#include "resfitD.C"
// FIT QUALITY under the adopted configuration, every theta_CM bin, both codes.
//   vz in [300,900] · TIE12 ON · set D12 + 6.68 · continuum p3* q / Ex^n at the corrected Ex zero
// Data with Poisson errors, the total model, and each component drawn separately, so the fit can
// be judged by eye and not only by a chi2 number.
double compAt(int b,double X,int which){            // 0 tot 1 cont 2 D12 3 6.68 4 gs+5/2
  RK::Bin& d=RK::B[b];
  double Xt=X-d.gsMu;
  double gs = d.gsAmp*shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
            + d.e5Amp*shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
  double bg = RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
  double d12=(d.in[0]?RK::LASTA[b][0]*d.T[0].Eval(X):0.0)
            +(d.in[1]?RK::LASTA[b][1]*d.T[1].Eval(X):0.0);
  double d34=(d.in[2]?RK::LASTA[b][2]*d.T[2].Eval(X):0.0)
            +(d.in[3]?RK::LASTA[b][3]*d.T[3].Eval(X):0.0);
  double e37=(d.in[6]?RK::LASTA[b][6]*d.T[6].Eval(X):0.0);
  switch(which){ case 1: return bg; case 2: return d12; case 3: return d34;
                 case 4: return gs; case 5: return e37; }
  return gs+bg+d12+d34+e37;
}
void fitq7(){
  setCM(10,5,60);
  RK::exLo=3.0; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true;            // 4.43/4.93 and 5.98/6.68 each tied, global ratio
  RK::EPP6[6]=3.70+RK::PPGS6;                 // the seventh component the residuals demanded
  bool keep[7]={true,true,true,true,false,false,true};
  for(int i=0;i<RK::NR;i++) RK::USE[i]=keep[i];
  Cfg b0=CFG[0], b1=CFG[1];
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",2000,1500); c->Divide(5,4,0.001,0.001);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
    rkRun(CFG[ic],NCM,false);
    printf("\n==== %s : fit quality, adopted configuration ====\n",CFG[ic].nm);
    printf("   %-10s %8s %9s %7s\n","theta_CM","N","chi2/ndf","ndf");
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      c->cd(1+ic*10+b); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone(Form("q%d%d",ic,b)); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(RK::exLo,RK::exHi);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=RK::exLo&&x<=RK::exHi) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(4.0*std::max(1.0,mx)); hd->SetMinimum(0.35);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.5); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ;E_{x} [MeV];counts / bin",
                        CFG[ic].nm,d.cmLo,d.cmHi));
      hd->GetXaxis()->SetTitleSize(0.055); hd->GetYaxis()->SetTitleSize(0.055);
      hd->GetXaxis()->SetLabelSize(0.05);  hd->GetYaxis()->SetLabelSize(0.05);
      hd->SetTitleSize(0.06);
      hd->Draw("E");
      // components
      struct K{int w; int col; int sty; int wid;};
      K ks[5]={{1,kGray+2,7,2},{5,kMagenta+2,1,2},{2,kGreen+3,1,2},{3,kOrange+7,1,2},{0,col[ic],1,3}};
      for(auto&k:ks){ auto* g=new TGraph();
        for(double x=RK::exLo;x<=RK::exHi+1e-9;x+=0.02) g->SetPoint(g->GetN(),x,compAt(b,x,k.w));
        g->SetLineColor(k.col); g->SetLineStyle(k.sty); g->SetLineWidth(k.wid); g->Draw("L"); }
      hd->Draw("E SAME");
      // chi2 on bins with >=5 counts
      double s=0; int nd=0;
      for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i),n=hd->GetBinContent(i);
        if(x<RK::exLo||x>RK::exHi||n<5) continue;
        double m=compAt(b,x,0); if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
      int ndf=std::max(1,nd-4);
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.062);
      tx->DrawLatex(0.56,0.86,Form("#chi^{2}/ndf %.2f",s/ndf));
      tx->SetTextSize(0.050); tx->SetTextColor(kGray+2);
      tx->DrawLatex(0.56,0.79,Form("N %.0f",d.N));
      printf("   %3.0f-%-6.0f %8.0f %9.2f %7d\n",d.cmLo,d.cmHi,d.N,s/ndf,ndf); }
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/fit_quality_7comp.png");
  printf("\n   wrote fit_quality_7comp.png\n");
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::TIE12=false;
  RK::PAIRS=false; for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6); BB_QUIET=true; binbuild();
}
