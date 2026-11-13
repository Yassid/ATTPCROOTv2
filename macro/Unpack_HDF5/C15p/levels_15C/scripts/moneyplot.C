#include "tiedfit.C"
// The result, on one panel: the two structures the data DOES determine, both codes, each
// normalised to its own maximum so the SHAPES can be compared directly.
void moneyplot(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.045,"t");
  double A[2][4][5], th[5];
  for(int b=0;b<5;b++) th[b]=TCM0+(b+0.5)*TCMSTEP;
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic];
    SF::B.clear(); SF::UCM.clear();
    for(int ib=0;ib<SF::NBIN;ib++) sfPrep(k,ib);
    for(int ib=0;ib<SF::NBIN;ib++) SF::UCM.push_back(std::cos(th[ib]*M_PI/180.0));
    SF::SMOOTH=false; SF::TIED=true;
    int npar=2+SF::NBIN*8+2;
    auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000); m->SetTolerance(0.1); m->SetPrintLevel(0);
    ROOT::Math::Functor fnc(&SF::fcn,npar); m->SetFunction(fnc);
    m->SetLimitedVariable(0,"bkg_w",1.0,0.02,0.05,3.0);
    m->SetLimitedVariable(1,"bkg_kappa",0.3,0.01,0.0,2.0);
    for(int b=0;b<SF::NBIN;b++){ SF::BinData&d=SF::B[b];
      double pk=d.h->GetBinContent(d.h->GetMaximumBin()); double lvl=0;int nl=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){double x=d.h->GetBinCenter(i);
        if(x>6&&x<8){lvl+=d.h->GetBinContent(i);nl++;}} lvl=nl?lvl/nl:1;
      int o=2+b*8;
      m->SetLimitedVariable(o+0,Form("gs_%d",b),pk,pk*1e-3,0,pk*3);
      m->SetLimitedVariable(o+1,Form("e52_%d",b),pk*0.05,pk*1e-4,0,pk);
      m->SetLimitedVariable(o+2,Form("cont_%d",b),1e-3,1e-5,0,1e3);
      m->SetLimitedVariable(o+3,Form("bkgB_%d",b),0.5*lvl,lvl*1e-3,0,lvl*20);
      for(int e=0;e<4;e++){ int a=SF::EMAP[e][0];
        double sc=d.in[a]?d.T[a].Eval((EPP[a]-PPGS)+d.gsMu):0; double seed=(sc>0)?0.2*lvl/sc:0;
        m->SetLimitedVariable(o+4+e,Form("E%d_%d",e+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*200+1e3);} }
    m->SetLimitedVariable(2+SF::NBIN*8+0,"r1",1.0,0.02,0.0,10.0);
    m->SetLimitedVariable(2+SF::NBIN*8+1,"r2",0.4,0.02,0.0,10.0);
    m->Minimize(); const double* p=m->X();
    for(int b=0;b<5;b++){ const double* q=p+2+b*8;
      for(int e=0;e<4;e++){ int a=SF::EMAP[e][0];
        double ph=SF::B[b].in[a]?q[4+e]*SF::B[b].T[a].Eval((EPP[a]-PPGS)+SF::B[b].gsMu):0;
        A[ic][e][b]=ph/std::sin(th[b]*M_PI/180.0); } }
    SF::TIED=false; delete m; }

  TCanvas* c=new TCanvas("cm","",1000,720);
  gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.09);
  auto* fr=gPad->DrawFrame(8,0.02,62,9.0);
  fr->SetTitle("^{15}C(d,d'): angular distributions of the two determined structures"
               ";#theta_{CM} (deg);d#sigma/d#Omega   (normalised to each maximum)");
  fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.25);
  // E1 green, E3 orange; ATTPCROOT filled circles + solid, Spyral open squares + dashed
  int col[2]={kGreen+2,kOrange+7}; int ei[2]={0,2};
  const char* enm[2]={"E1   4.43 + 4.93 MeV","E3   6.68 + 7.38 MeV"};
  // legend ABOVE the curves (frame extended to 9) so it cannot sit on a data point --
  // it was hiding the ATTPCROOT E3 point at 35 deg
  TLegend* lg=new TLegend(0.40,0.72,0.95,0.90);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.031);
  for(int j=0;j<2;j++){
    for(int ic=0;ic<2;ic++){
      double mx=0; for(int b=0;b<5;b++) mx=std::max(mx,A[ic][ei[j]][b]);
      if(mx<=0) continue;
      auto* g=new TGraph();
      for(int b=0;b<5;b++){ double v=A[ic][ei[j]][b]/mx; if(v>0.015) g->SetPoint(g->GetN(),th[b],v); }
      g->SetLineColor(col[j]); g->SetMarkerColor(col[j]); g->SetLineWidth(3);
      g->SetLineStyle(ic?2:1); g->SetMarkerStyle(ic?25:20); g->SetMarkerSize(1.8);
      g->Draw("LP");
      lg->AddEntry(g,Form("%s  %s",enm[j],ic?"(Spyral)":"(ATTPCROOT)"),"lp"); } }
  lg->Draw();
  TLatex t; t.SetNDC(true); t.SetTextSize(0.034);
  t.DrawLatex(0.17,0.30,"#color[418]{E1} peaks at the most forward angle;");
  t.DrawLatex(0.17,0.25,"#color[807]{E3} peaks one bin out, in BOTH codes.");
  t.DrawLatex(0.17,0.20,"Different shapes #Rightarrow different L transfer.");
  t.SetTextSize(0.029);
  t.DrawLatex(0.17,0.145,"5 angular points; the ordering is consistent, the shapes are not yet a measurement.");
  c->SaveAs("/home/yassid/Data/C15d_statefit/angdist_result.png");
  printf("wrote angdist_result.png\n");
}
