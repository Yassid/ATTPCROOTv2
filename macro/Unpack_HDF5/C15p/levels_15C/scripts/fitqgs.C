#include "resfitD.C"
// The same per-angle fits, drawn over the FULL range so the ground state and the 5/2+ are in view.
// The fit acts only above exLo (shaded edge); the g.s. and 5/2+ amplitudes and the response are
// measured on the full spectrum in rkPrep and then held fixed, so the curve below the window is a
// PREDICTION, not a fit -- which makes it a useful check in its own right.
//   adopted: vz [300,900] · chi2 loosened to 5.0 (ATTPCROOT) · both pairs tied, global ratio
double cAt(int b,double X,int which){   // 0 tot 1 cont 2 P1 3 P2 4 gs 5 5/2+ 6 singles
  RK::Bin& d=RK::B[b];
  double Xt=X-d.gsMu;
  double gs = d.gsAmp*shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
  double e5 = d.e5Amp*shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
  double bg = RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
  double p1=(d.in[0]?RK::LASTA[b][0]*d.T[0].Eval(X):0.0)+(d.in[1]?RK::LASTA[b][1]*d.T[1].Eval(X):0.0);
  double p2=(d.in[2]?RK::LASTA[b][2]*d.T[2].Eval(X):0.0)+(d.in[3]?RK::LASTA[b][3]*d.T[3].Eval(X):0.0);
  double sg=(d.in[4]?RK::LASTA[b][4]*d.T[4].Eval(X):0.0)+(d.in[5]?RK::LASTA[b][5]*d.T[5].Eval(X):0.0);
  switch(which){ case 1:return bg; case 2:return p1; case 3:return p2;
                 case 4:return gs; case 5:return e5; case 6:return sg; }
  return gs+e5+bg+p1+p2+sg;
}
void fitqgs(double exLo=3.5){
  setCM(10,5,60);
  RK::exLo=exLo; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::NPERBIN=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);        // the full six, as two tied pairs + 2 singles
  Cfg b0=CFG[0], b1=CFG[1];
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",2000,1500); c->Divide(5,4,0.001,0.001);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    rkRun(CFG[ic],NCM,false);
    printf("\n==== %s : fit quality, full range drawn, fit from %.1f ====\n",CFG[ic].nm,exLo);
    printf("   %-10s %9s %9s %9s\n","theta_CM","N","chi2/ndf","gs amp");
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      c->cd(1+ic*10+b); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone(Form("q%d%d",ic,b)); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(-2.0,11.0);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=-2&&x<=11) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(6.0*std::max(1.0,mx)); hd->SetMinimum(0.35);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.45); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ;E_{x} [MeV];counts / bin",
                        CFG[ic].nm,d.cmLo,d.cmHi));
      hd->GetXaxis()->SetTitleSize(0.055); hd->GetYaxis()->SetTitleSize(0.055);
      hd->GetXaxis()->SetLabelSize(0.05);  hd->GetYaxis()->SetLabelSize(0.05);
      hd->Draw("E");
      // shade the region the fit does NOT see
      { auto* bx=new TBox(-2.0,0.35,exLo,hd->GetMaximum()); bx->SetFillColorAlpha(kGray+1,0.22);
        bx->SetLineColor(0); bx->Draw(); hd->Draw("E SAME"); }
      struct K{int w;int col;int sty;int wid;};
      K ks[6]={{4,kViolet+1,1,2},{5,kCyan+2,1,2},{1,kGray+2,7,2},
               {2,kGreen+3,1,2},{3,kOrange+7,1,2},{0,col[ic],1,3}};
      for(auto&kk:ks){ auto* g=new TGraph();
        for(double x=-2.0;x<=11.0;x+=0.02){ double v=cAt(b,x,kk.w); if(v>0) g->SetPoint(g->GetN(),x,v); }
        if(g->GetN()<2) continue;
        g->SetLineColor(kk.col); g->SetLineStyle(kk.sty); g->SetLineWidth(kk.wid); g->Draw("L"); }
      hd->Draw("E SAME");
      double s=0; int nd=0;
      for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i),n=hd->GetBinContent(i);
        if(x<exLo||x>RK::exHi||n<5) continue;
        double m=cAt(b,x,0); if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
      int ndf=std::max(1,nd-4);
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.060);
      tx->DrawLatex(0.55,0.87,Form("#chi^{2}/ndf %.2f",s/ndf));
      tx->SetTextSize(0.048); tx->SetTextColor(kGray+2);
      tx->DrawLatex(0.55,0.80,Form("N %.0f",d.N));
      printf("   %3.0f-%-6.0f %9.0f %9.2f %9.0f\n",d.cmLo,d.cmHi,d.N,s/ndf,d.gsAmp); }
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/fit_quality_with_gs.png");
  printf("\n   wrote fit_quality_with_gs.png\n");
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
