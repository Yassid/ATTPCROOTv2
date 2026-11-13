#include "tiedfit.C"
// Angular distributions of the RESOLVABLE structures, both codes.
void tiedfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.05,"t");
  double A[2][4][5]; double th[5];
  for(int b=0;b<5;b++) th[b]=TCM0+(b+0.5)*TCMSTEP;
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic];
    SF::B.clear(); SF::UCM.clear();
    for(int ib=0;ib<SF::NBIN;ib++) sfPrep(k,ib);
    for(int ib=0;ib<SF::NBIN;ib++) SF::UCM.push_back(std::cos(th[ib]*M_PI/180.0));
    SF::SMOOTH=false; SF::TIED=true;
    int npar=2+SF::NBIN*8+2;
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000);
    m->SetTolerance(0.1); m->SetPrintLevel(0);
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
  TCanvas* c=new TCanvas("ct","",1500,640); c->Divide(2,1);
  const char* nm[4]={"E1  4.43+4.93","E2  5.98","E3  6.68+7.38","E4  8.28"};
  int col[4]={kGreen+2,kCyan+2,kOrange+7,kMagenta+2};
  for(int ic=0;ic<2;ic++){
    c->cd(ic+1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
    auto* fr=gPad->DrawFrame(8,0.3,62,400);
    fr->SetTitle(Form("%s  --  tied-doublet angular distributions"
                      ";#theta_{CM} (deg);amplitude / sin#theta   (#propto d#sigma/d#Omega)",CFG[ic].nm));
    fr->GetXaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetTitleSize(0.047);
    TLegend* lg=new TLegend(0.60,0.68,0.96,0.90); lg->SetBorderSize(0); lg->SetTextSize(0.040);
    for(int e=0;e<4;e++){ auto* g=new TGraph();
      for(int b=0;b<5;b++) if(A[ic][e][b]>0.3) g->SetPoint(g->GetN(),th[b],A[ic][e][b]);
      if(!g->GetN()) continue;
      g->SetLineColor(col[e]); g->SetMarkerColor(col[e]); g->SetLineWidth(3);
      g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->Draw("LP");
      lg->AddEntry(g,nm[e],"lp"); }
    lg->Draw();
    TLatex t; t.SetNDC(true); t.SetTextSize(0.037);
    t.DrawLatex(0.17,0.24,"E1 peaks at the most forward angle;");
    t.DrawLatex(0.17,0.18,"E3 peaks one bin out, in BOTH codes.");
  }
  c->SaveAs("/home/yassid/Data/C15d_statefit/tied_angdist.png");
  printf("wrote tied_angdist.png\n");
}
