#include "simfit.C"
// The angle-to-angle disagreement is in the DECOMPOSITION, not in the total.
// Left  : total inelastic strength per bin / sin(theta)  -- smooth, forward-peaked, both codes.
// Right : the six individual BW amplitudes -- erratic, alternating with zeros.
void angfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.05,"t");
  int col[6]={kGreen+2,kSpring-6,kCyan+2,kOrange+7,kRed-4,kMagenta+2};
  int mcol[2]={kAzure+2,kRed+1};
  TCanvas* c=new TCanvas("ca","",1650,680); c->Divide(2,1);
  std::vector<std::vector<double>> TOT(2), AMP[2];
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic];
    sfRun(k,false,false);
    SF::SMOOTH=false;
    int nper=SF::NPB, npar=2+SF::NBIN*nper;
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(1000000); m->SetMaxIterations(200000);
    m->SetTolerance(0.1); m->SetPrintLevel(0);
    ROOT::Math::Functor fnc(&SF::fcn,npar); m->SetFunction(fnc);
    m->SetLimitedVariable(0,"bkg_w",1.0,0.02,0.05,3.0);
    m->SetLimitedVariable(1,"bkg_kappa",0.3,0.01,0.0,2.0);
    for(int b=0;b<SF::NBIN;b++){ SF::BinData&d=SF::B[b]; if(!d.h) continue;
      double pk=d.h->GetBinContent(d.h->GetMaximumBin()); double lvl=0;int nl=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){double x=d.h->GetBinCenter(i);
        if(x>6&&x<8){lvl+=d.h->GetBinContent(i);nl++;}} lvl=nl?lvl/nl:1;
      int o=2+b*nper;
      m->SetLimitedVariable(o+0,Form("gs_%d",b),pk,pk*1e-3,0,pk*3);
      m->SetLimitedVariable(o+1,Form("e52_%d",b),pk*0.05,pk*1e-4,0,pk);
      m->SetLimitedVariable(o+2,Form("cont_%d",b),1e-3,1e-5,0,1e3);
      m->SetLimitedVariable(o+3,Form("bkgB_%d",b),0.5*lvl,lvl*1e-3,0,lvl*20);
      for(int i=0;i<6;i++){ double sc=d.in[i]?d.T[i].Eval((EPP[i]-PPGS)+d.gsMu):0;
        double seed=(sc>0)?0.2*lvl/sc:0;
        if(!d.in[i]){m->SetFixedVariable(o+4+i,Form("bw%d_%d",i+1,b),0);continue;}
        m->SetLimitedVariable(o+4+i,Form("bw%d_%d",i+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*200+1e3);} }
    m->Minimize(); const double* p=m->X();
    AMP[ic].assign(6,std::vector<double>());
    for(int b=0;b<SF::NBIN;b++){ SF::BinData&d=SF::B[b];
      double s=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i);
        if(x<3.5||x>9.0) continue;
        for(int j=0;j<6;j++) if(d.in[j]) s+=SF::ampOf(b,j,p)*d.T[j].Eval(x); }
      double th=(TCM0+(b+0.5)*TCMSTEP)*M_PI/180.0;
      TOT[ic].push_back(s/std::sin(th));
      for(int j=0;j<6;j++)
        AMP[ic][j].push_back(d.in[j]?SF::ampOf(b,j,p)*d.T[j].Eval((EPP[j]-PPGS)+d.gsMu):0); }
    delete m;
  }
  // ---- left: the total ----
  c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
  auto* fr=gPad->DrawFrame(5,200,65,20000);
  fr->SetTitle("TOTAL inelastic strength: smooth and consistent"
               ";#theta_{CM} (deg);#Sigma BW area / sin#theta   (#propto d#sigma/d#Omega)");
  fr->GetXaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetTitleSize(0.047);
  TLegend* lg=new TLegend(0.60,0.74,0.96,0.90); lg->SetBorderSize(0); lg->SetTextSize(0.042);
  for(int ic=0;ic<2;ic++){ auto* g=new TGraph();
    for(size_t b=0;b<TOT[ic].size();b++) g->SetPoint(b,TCM0+(b+0.5)*TCMSTEP,TOT[ic][b]);
    g->SetLineColor(mcol[ic]); g->SetMarkerColor(mcol[ic]); g->SetLineWidth(3);
    g->SetMarkerStyle(ic?21:20); g->SetMarkerSize(1.6); g->Draw("LP");
    lg->AddEntry(g,CFG[ic].nm,"lp"); }
  lg->Draw();
  TLatex t; t.SetNDC(true); t.SetTextSize(0.038);
  t.DrawLatex(0.17,0.24,"falls 7.3x (ATTPCROOT) and 4.1x (Spyral)");
  t.DrawLatex(0.17,0.18,"from 15 to 45 deg, same shape");
  // ---- right: the decomposition ----
  c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
  auto* f2=gPad->DrawFrame(5,0,65,110);
  f2->SetTitle("the DECOMPOSITION: not determined"
               ";#theta_{CM} (deg);BW peak height (counts / 100 keV)");
  f2->GetXaxis()->SetTitleSize(0.047); f2->GetYaxis()->SetTitleSize(0.047);
  for(int j=0;j<6;j++){ auto* g=new TGraph();
    for(size_t b=0;b<AMP[0][j].size();b++) g->SetPoint(b,TCM0+(b+0.5)*TCMSTEP,AMP[0][j][b]);
    g->SetLineColor(col[j]); g->SetMarkerColor(col[j]); g->SetLineWidth(3);
    g->SetMarkerStyle(20); g->SetMarkerSize(1.3); g->Draw("LP"); }
  TLatex t2; t2.SetNDC(true); t2.SetTextSize(0.040);
  for(int j=0;j<6;j++){ t2.SetTextColor(col[j]);
    t2.DrawLatex(0.63,0.88-j*0.05,Form("BW%d  %.2f MeV",j+1,EPP[j]-PPGS)); }
  t2.SetTextColor(kBlack); t2.SetTextSize(0.036);
  t2.DrawLatex(0.17,0.30,"ATTPCROOT. Amplitudes alternate with");
  t2.DrawLatex(0.17,0.25,"exact zeros and jump by factors of 3-7:");
  t2.DrawLatex(0.17,0.20,"six overlapping resonances trade freely");
  t2.DrawLatex(0.17,0.15,"inside each bin.");
  c->SaveAs("/home/yassid/Data/C15d_statefit/angular_decomposition.png");
  printf("wrote angular_decomposition.png\n");
}
