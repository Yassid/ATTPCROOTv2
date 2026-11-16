#include "resfitC.C"
// The continuum exponent against ONE variable: the low-vz cut, with vzHi fixed at 1000 for both
// codes so the only thing changing is where the window starts.  If n were a property of the
// physics it would sit flat on the s-wave halo value of 2.  It does not.
void nscanfig(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg b0=CFG[0], b1=CFG[1];
  CFG[0].vzHi=1000; CFG[1].vzHi=1000;
  double cuts[6]={0,100,200,300,400,500};
  TGraphErrors* g[2]; for(int i=0;i<2;i++) g[i]=new TGraphErrors();
  printf("   %-8s %12s %12s\n","VZLO","ATTPCROOT n","Spyral n");
  for(double v:cuts){
    VZLO=v; BB_QUIET=true; binbuild();
    double n[2],s[2];
    for(int ic=0;ic<2;ic++){
      RH::BKG=RH::HALON; double vf=rhRun(CFG[ic],NCM,false); n[ic]=RH::LASTP0;
      s[ic]=0; RH::SBFIX=true;
      for(double d=0.02;d<1.5;d+=0.02){ RH::SBVAL=n[ic]+d;
        if(rhRun(CFG[ic],NCM,false)-vf>=1.0){ s[ic]=d; break; } }
      RH::SBFIX=false;
      int p=g[ic]->GetN(); g[ic]->SetPoint(p,v,n[ic]); g[ic]->SetPointError(p,0,s[ic]); }
    printf("   %-8.0f %6.3f+-%.3f %6.3f+-%.3f\n",v,n[0],s[0],n[1],s[1]); }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; BB_QUIET=true; binbuild();
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",760,560);
  gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(-40,0.6,540,3.4);
  fr->SetTitle("continuum form factor  p_{3}*(M) q(M) / E_{x}^{n},  v_{z}^{max} = 1000 mm;"
               "low-v_{z} cut [mm];fitted exponent  n");
  { auto* l=new TLine(-40,2,540,2); l->SetLineColor(kGreen+3); l->SetLineWidth(3);
    l->SetLineStyle(2); l->Draw(); }
  int col[2]={kAzure+2,kRed+1}; int mk[2]={20,21};
  for(int ic=0;ic<2;ic++){ g[ic]->SetMarkerStyle(mk[ic]); g[ic]->SetMarkerSize(1.3);
    g[ic]->SetMarkerColor(col[ic]); g[ic]->SetLineColor(col[ic]); g[ic]->SetLineWidth(2);
    g[ic]->Draw("PL"); }
  auto* lg=new TLegend(0.17,0.73,0.62,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(g[0],"ATTPCROOT","pl"); lg->AddEntry(g[1],"Spyral","pl");
  lg->Draw();
  auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.033);
  tx->SetTextColor(kGreen+3); tx->DrawLatex(0.63,0.455,"s-wave halo: n = 2");
  tx->SetTextColor(kGray+3); tx->SetTextSize(0.030);
  tx->DrawLatex(0.17,0.185,"#splitline{error bars are STATISTICAL (#pm0.1).  The swing with the}"
                           "{acceptance window is ~10x larger: n is not measured.}");
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/continuum_exponent_unstable.png");
  printf("   wrote continuum_exponent_unstable.png\n");
}
