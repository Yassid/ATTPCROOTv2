#include "binacc.h"
#include "binbuild.C"
// The vz story in one figure: where the events are, where the adopted windows were, and where the
// Ex calibration actually holds.
void vzfig(){
  gStyle->SetOptStat(0);
  const int NS=10; double lo=0,hi=1000;
  TH1D* hv[2]; TGraphErrors* gc[2]; TGraphErrors* gf[2];
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic]; Pars A=parsOf(k); A.vzLo=-1e9; A.vzHi=1e9;
    hv[ic]=new TH1D(Form("hv%d",ic),"",110,-100,1100); hv[ic]->SetDirectory(nullptr);
    TH1D* h[NS]; for(int s=0;s<NS;s++){ h[s]=new TH1D(Form("q%d%d",ic,s),"",160,-2,2); h[s]->SetDirectory(nullptr); }
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||icv<k.icLo||icv>k.icHi) continue;
      hv[ic]->Fill(vz);
      if(vz<lo||vz>=hi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double e=exOf(A,ebAt(A,vz,true),thCorr(A,thf,ke),ke); if(!std::isfinite(e)) continue;
      int idx=(int)((vz-lo)/((hi-lo)/NS)); if(idx<0||idx>=NS) continue; h[idx]->Fill(e); }
    f->Close();
    gc[ic]=new TGraphErrors(); gf[ic]=new TGraphErrors();
    for(int s=0;s<NS;s++){ if(h[s]->GetEntries()<200) continue;
      TF1 g(Form("gz%d%d",ic,s),"gaus",-1.0,1.0);
      g.SetParameters(h[s]->GetBinContent(h[s]->GetMaximumBin()),0,0.25);
      h[s]->Fit(&g,"RQ0N");
      double x=lo+(s+0.5)*(hi-lo)/NS;
      int p=gc[ic]->GetN();
      gc[ic]->SetPoint(p,x,g.GetParameter(1)); gc[ic]->SetPointError(p,50,g.GetParError(1));
      gf[ic]->SetPoint(p,x,2.355*g.GetParameter(2)); gf[ic]->SetPointError(p,50,2.355*g.GetParError(2)); }
    for(auto*g:{gc[ic],gf[ic]}){ g->SetMarkerStyle(20); g->SetMarkerSize(1.1);
      g->SetMarkerColor(col[ic]); g->SetLineColor(col[ic]); g->SetLineWidth(2); }
  }
  TCanvas* c=new TCanvas("c","c",1500,470); c->Divide(3,1);
  // (1) vz distributions with the adopted windows
  c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
  hv[0]->SetTitle(";reconstructed vertex v_{z} [mm];gated (d,d') events / 10 mm");
  hv[0]->SetLineColor(col[0]); hv[0]->SetLineWidth(2);
  hv[1]->SetLineColor(col[1]); hv[1]->SetLineWidth(2);
  hv[0]->SetMaximum(1.35*std::max(hv[0]->GetBinContent(hv[0]->GetMaximumBin()),
                                  hv[1]->GetBinContent(hv[1]->GetMaximumBin())));
  hv[0]->Draw("HIST"); hv[1]->Draw("HIST SAME");
  { auto* b=new TBox(0,0,300,hv[0]->GetMaximum()); b->SetFillColorAlpha(kGray+1,0.30);
    b->SetLineColor(0); b->Draw(); hv[0]->Draw("HIST SAME"); hv[1]->Draw("HIST SAME"); }
  { auto* l=new TLine(400,0,400,hv[0]->GetMaximum()*0.78); l->SetLineColor(col[1]);
    l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    auto* tx=new TLatex(); tx->SetTextSize(0.034); tx->SetTextColor(col[1]);
    tx->DrawLatex(410,hv[0]->GetMaximum()*0.80,"Spyral v_{z}^{max}=400"); }
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.034);
    tx->SetTextColor(col[0]); tx->DrawLatex(0.55,0.86,"ATTPCROOT  (0-1100)");
    tx->SetTextColor(col[1]); tx->DrawLatex(0.55,0.80,"Spyral  (0-400)");
    tx->SetTextColor(kGray+3); tx->SetTextSize(0.030);
    tx->DrawLatex(0.20,0.60,"#splitline{shaded: v_{z} < 300,}{calibration breaks}"); }
  // (2) centroid
  c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(0,-0.45,1000,0.45);
  fr->SetTitle("ground-state centroid;v_{z} [mm];E_{x} centroid [MeV]");
  { auto* b=new TBox(0,-0.45,300,0.45); b->SetFillColorAlpha(kGray+1,0.30); b->SetLineColor(0); b->Draw(); }
  for(int ic=0;ic<2;ic++) gc[ic]->Draw("P");
  { auto* l=new TLine(0,0,1000,0); l->SetLineStyle(2); l->Draw(); }
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.034); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.40,0.88,"walks 0.212 / 0.123 MeV below 300"); }
  // (3) FWHM
  c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr2=gPad->DrawFrame(0,0.35,1000,1.30);
  fr2->SetTitle("ground-state resolution;v_{z} [mm];E_{x} FWHM [MeV]");
  { auto* b=new TBox(0,0.35,300,1.30); b->SetFillColorAlpha(kGray+1,0.30); b->SetLineColor(0); b->Draw(); }
  for(int ic=0;ic<2;ic++) gf[ic]->Draw("P");
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.034); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.33,0.88,"peak broadens 10-30 % below 300"); }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/calibration_vs_vz.png");
  printf("   wrote calibration_vs_vz.png\n");
}
