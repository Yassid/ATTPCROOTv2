#include "resfitC.C"
// The continuum, isolated.  Points = data minus everything else the fit knows about (the
// resonances and the g.s./5/2+ tails).  Curves = the candidate continuum models, each at its own
// best fit.  Uniform three-body phase space is the WRONG SHAPE; multiplying it by the s-wave halo
// momentum distribution 1/Ex^2 -- no free parameter -- fixes it.
TGraph* curveOf(int b,int form,double par,double sc){
  auto* g=new TGraph(); RH::Bin& d=RH::B[b];
  int sv=RH::BKG; RH::BKG=form;
  for(double x=RH::exLo;x<=RH::exHi+1e-6;x+=0.02){
    double a=d.ga?std::max(0.0,d.ga->Eval(x)):1.0;
    g->SetPoint(g->GetN(),x,sc*RH::bshape(x,par)*a); }
  RH::BKG=sv; return g;
}
void bkmodfig(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  VZLO=200; BB_QUIET=true; binbuild();   // remove the vz 90-180 mm instrumental excess
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",1150,520); c->Divide(2,1);
  const int IB=0;   // theta_CM 10-20
  TLegend* lg=nullptr;
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic];
    RH::BKG=RH::HALON; rhRun(k,NCM,false); double nhat=RH::LASTP0;
    RH::Bin& d=RH::B[IB];
    // points: data minus resonances and the g.s./5/2+ tails  (= the fitted continuum + noise)
    auto* gd=new TGraphErrors();
    for(int i=1;i<=d.h->GetNbinsX();i++){
      double x=d.h->GetBinCenter(i); if(x<RH::exLo||x>RH::exHi) continue;
      double n=d.h->GetBinContent(i);
      double other=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                  +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
      for(int r=0;r<RH::NR;r++) if(RH::USE[r]&&d.in[r]) other+=RH::LASTA[IB][r]*d.T[r].Eval(x);
      int p=gd->GetN(); gd->SetPoint(p,x,n-other); gd->SetPointError(p,0,std::sqrt(std::max(1.0,n))); }
    // each model at its own best fit, same bin
    struct M{int form; const char* nm; int col; int sty;};
    M ms[3]={{RH::PS,  "uniform 3-body  p_{3}*(M) q(M)",kRed+1,  2},
             {RH::HALO,"s-wave halo  p_{3}* q / E_{x}^{2}   (no free parameter)",kAzure+2,1},
             {RH::EXP, "empirical  exp(-#kappa E_{x})",kGray+2,7}};
    std::vector<TGraph*> cv;
    for(auto&m:ms){ RH::BKG=m.form; rhRun(k,NCM,false);
      double sc=0; { const int o=IB; sc=RH::LASTBG[o]; }
      auto* g=curveOf(IB,m.form,RH::LASTP0,sc);
      g->SetLineColor(m.col); g->SetLineWidth(3); g->SetLineStyle(m.sty); cv.push_back(g); }
    c->cd(ic+1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    gd->SetMarkerStyle(20); gd->SetMarkerSize(0.7);
    gd->SetTitle(Form("%s   #theta_{CM} 10-15#circ,  v_{z} > 200 mm;E_{x} [MeV];counts / bin  (data #minus resonances)",k.nm));
    gd->GetXaxis()->SetLimits(RH::exLo,RH::exHi);
    double mx=0; for(int i=0;i<gd->GetN();i++) mx=std::max(mx,gd->GetY()[i]);
    gd->SetMinimum(0); gd->SetMaximum(1.45*mx);
    gd->Draw("AP");
    for(auto*g:cv) g->Draw("L");
    gd->Draw("P");
    if(!lg){ lg=new TLegend(0.30,0.705,0.97,0.895); lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->SetTextSize(0.0285);
      for(size_t i=0;i<cv.size();i++) lg->AddEntry(cv[i],ms[i].nm,"l"); }
    lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.032);
    tx->DrawLatex(0.17,0.655,Form("p_{3}* q / E_{x}^{n},  n free:  n = %.2f",nhat));
    tx->DrawLatex(0.17,0.605,"s-wave halo predicts n = 2");
  }
  VZLO=0; BB_QUIET=true; binbuild();
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/continuum_model.png");
  printf("   wrote continuum_model.png\n");
}
