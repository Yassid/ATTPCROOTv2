// The acceptance cliff, measured on the existing Geant4 + reconstruction simulation.
// Flat (0.98-1.00) for vz 0-900, then an Ex-DEPENDENT fall in the last 100 mm.
// ATTPCROOT's adopted window (0-1100) includes that cliff; Spyral's (0-400) sits instead in the
// region where the Ex calibration has broken down.  The good window is vz in [300,900].
void accfig(){
  gStyle->SetOptStat(0);
  TFile f("/home/yassid/Data/C15dd_sim/acceptance3d_C15dd.root");
  std::vector<double> ex; std::vector<TString> tg;
  TIter it(f.GetListOfKeys()); TKey* k;
  while((k=(TKey*)it())){ TString n=k->GetName();
    if(!n.BeginsWith("eff_ex")) continue;
    tg.push_back(n); ex.push_back(TString(n(6,n.Length()-6)).Atof()); }
  std::vector<size_t> ord(ex.size()); for(size_t i=0;i<ord.size();i++) ord[i]=i;
  std::sort(ord.begin(),ord.end(),[&](size_t a,size_t b){return ex[a]<ex[b];});
  TCanvas* c=new TCanvas("c","c",1180,520); c->Divide(2,1);
  // (1) efficiency vs vz, one curve per Ex
  c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(0,0.25,1000,1.12);
  fr->SetTitle("PRA efficiency, simulation;true vertex v_{z} [mm];efficiency");
  auto* lg=new TLegend(0.17,0.17,0.46,0.52); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->SetTextSize(0.030); lg->SetHeader("E_{x} [MeV]");
  int nc=ord.size(); int ip=0;
  for(size_t j:ord){ auto* h=(TH2D*)f.Get(tg[j]); if(!h) continue;
    auto* hg=(TH2D*)f.Get("g_"+tg[j](4,tg[j].Length()-4));
    auto* hn=(TH2D*)f.Get("n_"+tg[j](4,tg[j].Length()-4));
    auto* g=new TGraph();
    for(int v=1;v<=h->GetNbinsY();v++){
      double G=hg?hg->Integral(1,hg->GetNbinsX(),v,v):0, N=hn?hn->Integral(1,hn->GetNbinsX(),v,v):0;
      if(G<=0) continue; g->SetPoint(g->GetN(),h->GetYaxis()->GetBinCenter(v),N/G); }
    int col=TColor::GetColorPalette((int)(254.0*ip/std::max(1,nc-1)));
    g->SetLineColor(col); g->SetLineWidth(3); g->SetMarkerColor(col);
    g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->Draw("PL");
    if(ip%2==0||ip==nc-1) lg->AddEntry(g,Form("%.2f",ex[j]),"l");
    ip++; }
  lg->Draw();
  { auto* b=new TBox(900,0.25,1000,1.12); b->SetFillColorAlpha(kRed,0.13); b->SetLineColor(0); b->Draw();
    auto* b2=new TBox(0,0.25,300,1.12); b2->SetFillColorAlpha(kGray+1,0.28); b2->SetLineColor(0); b2->Draw(); }
  gPad->RedrawAxis();
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.029);
    tx->SetTextColor(kRed+1);  tx->DrawLatex(0.60,0.30,"#splitline{acceptance}{cliff}");
    tx->SetTextColor(kGray+3); tx->DrawLatex(0.17,0.86,"#splitline{E_{x} calibration}{broken}"); }
  // (2) the cliff depth vs Ex -- the thing that distorts a spectrum
  c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* g2=new TGraph();
  for(size_t j:ord){ auto* hg=(TH2D*)f.Get("g_"+tg[j](4,tg[j].Length()-4));
    auto* hn=(TH2D*)f.Get("n_"+tg[j](4,tg[j].Length()-4)); if(!hg||!hn) continue;
    int v=hg->GetNbinsY();
    double G=hg->Integral(1,hg->GetNbinsX(),v,v), N=hn->Integral(1,hn->GetNbinsX(),v,v);
    if(G>0) g2->SetPoint(g2->GetN(),ex[j],N/G); }
  g2->SetMarkerStyle(20); g2->SetMarkerSize(1.4); g2->SetMarkerColor(kRed+1);
  g2->SetLineColor(kRed+1); g2->SetLineWidth(3);
  g2->SetTitle("efficiency in the last slice, v_{z} 900-1000 mm;E_{x} [MeV];efficiency");
  g2->GetXaxis()->SetLimits(-0.5,12.5); g2->SetMinimum(0.2); g2->SetMaximum(1.05);
  g2->Draw("APL");
  { auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.031); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.30,0.80,"#splitline{strongly E_{x}-DEPENDENT:}{0.86 at the g.s. #rightarrow 0.35 at 12 MeV}");
    tx->SetTextSize(0.028);
    tx->DrawLatex(0.17,0.22,"#splitline{ATTPCROOT's adopted window (0-1100) includes this,}{so its high-E_{x} continuum was suppressed uncorrected}"); }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/acceptance_cliff.png");
  printf("   wrote acceptance_cliff.png\n");
}
