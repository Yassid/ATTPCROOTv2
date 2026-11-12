#include "binfit.C"
// The resonance strengths against the choice of low-energy background model.
void bkgfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.05,"t");
  const char* MOD[7]={"none","erf @ S_{n}","erf @ 0.0","erf @ 1.0","erf @ 2.0","erf @ 2.5","long tail"};
  TCanvas* c=new TCanvas("cs","",1600,680); c->Divide(2,1);
  int col[NB]={kGreen+2,kSpring-6,kCyan+2,kOrange+7,kRed-4,kMagenta+2};
  for(int ic=0;ic<2;ic++){
    Cfg&k=CFG[ic];
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false; USE_LONGTAIL=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    std::vector<std::array<double,NB>> V; std::vector<int> ok;
    auto grab=[&](){ BR o=fitBin(k,0,-2.0,12.0,false,true);
      std::array<double,NB> a; for(int i=0;i<NB;i++) a[i]=o.amp[i];
      bool good=(o.bgK<0.45&&o.chi2ndf<6)||(!USE_LEBKG);
      V.push_back(a); ok.push_back(good?1:0); };
    USE_LEBKG=false; USE_LONGTAIL=false; FREE_TURNON=false; grab();
    USE_LEBKG=true;  grab();
    FREE_TURNON=true;
    for(double t : {0.0,1.0,2.0,2.5}){ TURNON=t; grab(); }
    FREE_TURNON=false; USE_LEBKG=false; USE_LONGTAIL=true; grab();
    USE_LONGTAIL=false;

    c->cd(ic+1); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.20); gPad->SetGridy();
    auto* fr=gPad->DrawFrame(-0.5,0,6.5,200);
    fr->SetTitle(Form("%s  #theta_{CM} 10-20#circ: resonance strength vs background model"
                      ";;BW peak height (counts / 100 keV)",k.nm));
    fr->GetYaxis()->SetTitleSize(0.045);
    for(int i=0;i<7;i++) fr->GetXaxis()->SetBinLabel(
        fr->GetXaxis()->FindBin(i),MOD[i]);
    fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetLabelSize(0.055);
    for(int i=0;i<NB;i++){
      auto* g=new TGraph();
      for(size_t j=0;j<V.size();j++){ if(!ok[j]) continue; g->SetPoint(g->GetN(),(double)j,V[j][i]); }
      g->SetLineColor(col[i]); g->SetMarkerColor(col[i]); g->SetLineWidth(3);
      g->SetMarkerStyle(20); g->SetMarkerSize(1.3); g->Draw("LP"); }
    TLatex t; t.SetNDC(true); t.SetTextSize(0.040);
    for(int i=0;i<NB;i++){ t.SetTextColor(col[i]);
      t.DrawLatex(0.66,0.88-i*0.05,Form("BW%d  %.2f MeV",i+1,EPP[i]-PPGS)); }
    t.SetTextColor(kBlack); t.SetTextSize(0.036);
    t.DrawLatex(0.16,0.90,"BW4 is stable; BW1-BW3 are not");
  }
  c->SaveAs("/home/yassid/Data/C15d_statefit/bkg_systematic.png");
  printf("wrote bkg_systematic.png\n");
}
