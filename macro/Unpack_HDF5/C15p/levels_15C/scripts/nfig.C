// The continuum exponent against the vz window.  Every point here was fitted in this session
// (scripts bkvzlo.C / vzgood.C / vzgood2.C / spyvz.C); error bars are STATISTICAL.
// If n were a property of the physics it would sit on the s-wave halo value of 2 regardless of
// which part of the chamber we accept.  It does not.
void nfig(){
  gStyle->SetOptStat(0);
  // label, ATTPCROOT n +- s, Spyral n +- s   (-1 = not measured in that configuration)
  struct P{const char* lab; double a,ea,s,es;};
  P pts[6]={
    {"adopted#splitline{}{(0-1100 / 0-400)}", 1.832,0.096, 1.220,0.077},
    {"0-900",                                 1.815,0.097, 2.147,0.059},
    {"0-1000",                                   -1,0,     2.147,0.059},
    {"200-1100 / 200-400",                    2.034,0.120, 1.867,0.120},
    {"300-1000",                              2.292,0.127, 2.879,0.087},
    {"300-900",                               2.262,0.129, 2.895,0.090}};
  const int N=6;
  TCanvas* c=new TCanvas("c","c",900,560);
  gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.22); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(-0.6,0.7,N-0.4,3.5);
  fr->SetTitle("continuum form factor  p_{3}*(M) q(M) / E_{x}^{n};;fitted exponent  n");
  fr->GetXaxis()->SetNdivisions(0);
  { auto* l=new TLine(-0.6,2,N-0.4,2); l->SetLineColor(kGreen+3); l->SetLineWidth(3);
    l->SetLineStyle(2); l->Draw(); }
  auto* ga=new TGraphErrors(); auto* gs=new TGraphErrors();
  for(int i=0;i<N;i++){
    if(pts[i].a>0){ int p=ga->GetN(); ga->SetPoint(p,i-0.08,pts[i].a); ga->SetPointError(p,0,pts[i].ea); }
    if(pts[i].s>0){ int p=gs->GetN(); gs->SetPoint(p,i+0.08,pts[i].s); gs->SetPointError(p,0,pts[i].es); }
    auto* tx=new TLatex(i,0.52,pts[i].lab); tx->SetTextAlign(22); tx->SetTextSize(0.028);
    tx->SetTextAngle(0); tx->Draw(); }
  ga->SetMarkerStyle(20); ga->SetMarkerSize(1.4); ga->SetMarkerColor(kAzure+2);
  ga->SetLineColor(kAzure+2); ga->SetLineWidth(3); ga->Draw("P");
  gs->SetMarkerStyle(21); gs->SetMarkerSize(1.4); gs->SetMarkerColor(kRed+1);
  gs->SetLineColor(kRed+1); gs->SetLineWidth(3); gs->Draw("P");
  auto* lg=new TLegend(0.16,0.76,0.44,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(ga,"ATTPCROOT","pl"); lg->AddEntry(gs,"Spyral","pl"); lg->Draw();
  auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.031);
  tx->SetTextColor(kGreen+3); tx->DrawLatex(0.70,0.42,"s-wave halo: n = 2");
  tx->SetTextColor(kGray+3); tx->SetTextSize(0.028);
  tx->DrawLatex(0.30,0.885,"error bars are STATISTICAL (#pm0.1); the swing with the v_{z} window is ~10#times larger");
  tx->SetTextSize(0.026);
  tx->DrawLatex(0.40,0.055,"accepted v_{z} window [mm]");
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/continuum_exponent_unstable.png");
  printf("   wrote continuum_exponent_unstable.png\n");
}
