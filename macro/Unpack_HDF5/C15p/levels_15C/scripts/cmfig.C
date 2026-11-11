#include "binbuild.C"
// The states stay at fixed Ex across theta_CM bins -- the decisive difference from theta_lab,
// where the Jacobian pile-up made every slice peak at its own kinematic edge.
void cmfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.07,"t");
  TFile f("dd_cmbins.root");
  TCanvas* c=new TCanvas("cc","",1750,900); c->Divide(5,2,0.001,0.001);
  int col[2]={kAzure+2,kRed+1};
  for(int r=0;r<2;r++){ const char* nm=CFG[r].nm;
    for(int i=0;i<5;i++){
      c->cd(r*5+i+1); gPad->SetLogy(); gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.02);
      gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.13);
      TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",nm,i)); if(!h) continue;
      h->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ",nm,TCM0+i*TCMSTEP,TCM0+(i+1)*TCMSTEP));
      h->GetXaxis()->SetRangeUser(-2,14); h->SetMinimum(0.5);
      h->SetLineColor(col[r]); h->SetLineWidth(2);
      h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
      h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitle("E_{x} (MeV)"); h->GetYaxis()->SetTitle("counts / 100 keV");
      h->Draw("hist");
      for(double e : {4.6,6.9}){ TLine* l=new TLine(e,0.5,e,h->GetMaximum());
        l->SetLineColor(e<5?kGreen+2:kOrange+7); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
      TLatex t; t.SetNDC(true); t.SetTextSize(0.055);
      t.DrawLatex(0.55,0.86,Form("N = %.0f",h->GetEntries()));
    }
  }
  c->SaveAs("/home/yassid/Data/C15d_statefit/dd_cmbins.png");
  printf("wrote dd_cmbins.png\n");
}
