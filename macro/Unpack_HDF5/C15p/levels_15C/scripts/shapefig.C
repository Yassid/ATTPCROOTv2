// generated from allang.log -- the overnight sweep
void shapefig(){
  gStyle->SetOptStat(0);
  const int N=10;
  double X[N]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  double aA[N]={1.166,1.993,0.932,0,0,0.5,0.75,0,0,0};
  double kA[N]={0.636,0.958,0.511,0.183,0.014,0,0,0,0,0};
  double dA[N]={632.8,899.8,492.2,2.6,-565.5,19,139.7,2.7,-42,-29.1};
  double cA[N]={1.43,1.76,2.46,2.24,14.92,1.15,0.78,1.1,1.62,1.15};
  double aS[N]={1.097,2.077,0.744,0.75,0.5,0.158,0.575,0,0.25,0};
  double kS[N]={0.686,1.004,0.506,0.305,0,0,0,0,8.969,1.748};
  double dS[N]={158.5,575.3,-2.2,-29.4,24.2,11.4,21.9,-14.1,222427,-268};
  double cS[N]={2.7,6.43,3.59,1.91,1.29,1.44,1.22,1.39,6.6,4.76};

  TCanvas* c=new TCanvas("c","c",1500,520); c->Divide(3,1);
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  // (1) the two shape parameters
  c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(8,-0.2,62,2.6);
  fr->SetTitle("background shape vs angle;#theta_{CM} [deg];a  (filled)   and   k  (open)");
  auto mk1=[&](double*Y,int ic,int st){ auto* g=new TGraph(N,X,Y);
    g->SetMarkerStyle(st); g->SetMarkerSize(1.4); g->SetMarkerColor(col[ic]);
    g->SetLineColor(col[ic]); g->SetLineWidth(2); g->Draw("PL"); return g; };
  auto* g1=mk1(aA,0,20); mk1(kA,0,24); auto* g2=mk1(aS,1,21); mk1(kS,1,25);
  auto* lg=new TLegend(0.52,0.72,0.92,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(g1,"ATTPCROOT","pl"); lg->AddEntry(g2,"Spyral","pl"); lg->Draw();
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.030); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.18,0.24,"#splitline{both #rightarrow 0 backward:}{PURE PHASE SPACE}"); }
  // (2) does the tailed peak beat a Gaussian?
  c->cd(2); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr2=gPad->DrawFrame(8,-700,62,1100);
  fr2->SetTitle("peak shape: tailed vs pure Gaussian;#theta_{CM} [deg];#Delta(-2lnL)   positive = tails BETTER");
  for(int ic=0;ic<2;ic++){ auto* g=new TGraph(N,X,ic?dS:dA);
    g->SetMarkerStyle(mk[ic]); g->SetMarkerSize(1.4); g->SetMarkerColor(col[ic]);
    g->SetLineColor(col[ic]); g->SetLineWidth(2); g->Draw("PL"); }
  { auto* l=new TLine(8,0,62,0); l->SetLineStyle(2); l->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.029); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.30,0.80,"#splitline{tails matter only where the}{g.s. peak is large (forward)}"); }
  // (3) fit quality
  c->cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr3=gPad->DrawFrame(8,0,62,7.5);
  fr3->SetTitle("fit quality, tailed peak, free background;#theta_{CM} [deg];#chi^{2}/ndf");
  for(int ic=0;ic<2;ic++){ auto* g=new TGraph(N,X,ic?cS:cA);
    g->SetMarkerStyle(mk[ic]); g->SetMarkerSize(1.4); g->SetMarkerColor(col[ic]);
    g->SetLineColor(col[ic]); g->SetLineWidth(2); g->Draw("PL"); }
  { auto* l=new TLine(8,1,62,1); l->SetLineStyle(2); l->SetLineColor(kGreen+3); l->Draw(); }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/shape_vs_angle.png");
  printf("   wrote shape_vs_angle.png\n");
}
