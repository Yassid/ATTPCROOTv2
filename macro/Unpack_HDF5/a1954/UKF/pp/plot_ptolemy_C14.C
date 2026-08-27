// Cross-checks of the a1954 14C FRESCO calculations against PtolemyCpp.
// Three separate figures, each panelled.
TGraph* rd(TString f,double sc=1){
  auto*g=new TGraph(); std::ifstream in(f.Data()); double a,b;
  while(in>>a>>b) if(b>0) g->SetPoint(g->GetN(),a,sc*b);
  if(!g->GetN()) printf("\033[1;31mempty %s\033[0m\n",f.Data());
  return g;
}
void plot_ptolemy_C14(){
  gStyle->SetOptStat(0); gStyle->SetTitleSize(0.055,"t");
  TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  TString P=here+"/../ptolemy/dat/";
  TString F=here+"/../fresco/outputs/";
  TString O="/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";

  // ---------- 1. elastic: the two codes on identical KD03 input ----------
  auto*c1=new TCanvas("c1","",1100,520); c1->Divide(2,1);
  c1->cd(1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
  auto*gf=rd(F+"p14C_el_161_dsdo.dat"); auto*gp=rd(P+"el_kd03_frescoparams.dat");
  gf->SetTitle("p+^{14}C elastic, identical KD03 input;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
  gf->SetLineColor(kBlue+1); gf->SetLineWidth(3); gf->Draw("AL");
  gp->SetLineColor(kRed); gp->SetLineWidth(2); gp->SetLineStyle(2); gp->Draw("L same");
  auto*l1=new TLegend(0.45,0.72,0.88,0.88); l1->SetBorderSize(0); l1->SetFillStyle(0);
  l1->AddEntry(gf,"FRESCO","l"); l1->AddEntry(gp,"PtolemyCpp","l"); l1->Draw();
  c1->cd(2); gPad->SetGridx(); gPad->SetGridy();
  auto*gr=new TGraph();
  for(int i=0;i<gp->GetN();i++){ double x=gp->GetX()[i]; if(x<5) continue;
    double f=gf->Eval(x); if(f>0) gr->SetPoint(gr->GetN(),x,gp->GetY()[i]/f); }
  gr->SetTitle("ratio Ptolemy / FRESCO;#theta_{cm} [deg];ratio");
  gr->SetLineColor(kBlack); gr->SetLineWidth(2);
  gr->SetMinimum(0.995); gr->SetMaximum(1.005); gr->Draw("AL");
  auto*u=new TLine(5,1,180,1); u->SetLineColor(kRed); u->SetLineStyle(2); u->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.042);
  t.DrawLatex(0.16,0.83,"agreement 0.05% -- Ptolemy prints 4 s.f.");
  c1->SaveAs(O+"01_elastic_fresco_vs_ptolemy.png");

  // ---------- 2. optical-potential sweep: does the dip move? ----------
  const int NP=5; const char* pk[NP]={"K","V","G","P","M"};
  const char* pn[NP]={"KD03 (Koning-Delaroche)","CH89 (Varner)","Becchetti-Greenlees","Perey","Menet"};
  int col[NP]={kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1};
  auto*c2=new TCanvas("c2","",1200,900); c2->Divide(2,2);
  // (a) all five, full range
  c2->cd(1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
  auto*lg2=new TLegend(0.34,0.62,0.93,0.89); lg2->SetBorderSize(0); lg2->SetFillStyle(0);
  lg2->SetTextSize(0.038);
  for(int i=0;i<NP;i++){ auto*g=rd(P+TString::Format("el_omp_%s.dat",pk[i]));
    g->SetLineColor(col[i]); g->SetLineWidth(2);
    if(i==0){ g->SetTitle("p+^{14}C elastic, five global proton OMPs;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
              g->GetXaxis()->SetLimits(0,180); g->Draw("AL"); } else g->Draw("L same");
    lg2->AddEntry(g,pn[i],"l"); }
  lg2->Draw();
  // (b) zoom on the dip
  c2->cd(2); gPad->SetGridx(); gPad->SetGridy();
  for(int i=0;i<NP;i++){ auto*g=rd(P+TString::Format("el_omp_%s.dat",pk[i]));
    g->SetLineColor(col[i]); g->SetLineWidth(3);
    if(i==0){ g->SetTitle("the diffraction minimum;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
              g->GetXaxis()->SetLimits(50,75); g->SetMinimum(0); g->SetMaximum(12); g->Draw("AL"); }
    else g->Draw("L same"); }
  TLatex t2; t2.SetNDC(); t2.SetTextSize(0.042);
  t2.DrawLatex(0.15,0.85,"dip walks 58#circ #rightarrow 63#circ");
  t2.DrawLatex(0.15,0.79,"depth varies #times2.6");
  // (c) the luminosity window
  c2->cd(3); gPad->SetGridx(); gPad->SetGridy();
  for(int i=0;i<NP;i++){ auto*g=rd(P+TString::Format("el_omp_%s.dat",pk[i]));
    g->SetLineColor(col[i]); g->SetLineWidth(3);
    if(i==0){ g->SetTitle("luminosity window (used for the normalisation);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
              g->GetXaxis()->SetLimits(70,140); g->SetMinimum(0); g->SetMaximum(22); g->Draw("AL"); }
    else g->Draw("L same"); }
  t2.DrawLatex(0.42,0.85,"KD03 sits high;");
  t2.DrawLatex(0.42,0.79,"the other four agree to #pm4%");
  // (d) ratio to KD03
  c2->cd(4); gPad->SetGridx(); gPad->SetGridy();
  auto*gK=rd(P+"el_omp_K.dat");
  for(int i=1;i<NP;i++){ auto*g=rd(P+TString::Format("el_omp_%s.dat",pk[i]));
    auto*q=new TGraph();
    for(int j=0;j<g->GetN();j++){ double x=g->GetX()[j]; if(x<10||x>150) continue;
      double k=gK->Eval(x); if(k>0) q->SetPoint(q->GetN(),x,g->GetY()[j]/k); }
    q->SetLineColor(col[i]); q->SetLineWidth(2);
    if(i==1){ q->SetTitle("ratio to KD03;#theta_{cm} [deg];#sigma / #sigma_{KD03}");
              q->SetMinimum(0); q->SetMaximum(3); q->GetXaxis()->SetLimits(10,150); q->Draw("AL"); }
    else q->Draw("L same"); }
  auto*u2=new TLine(10,1,150,1); u2->SetLineColor(kBlack); u2->SetLineStyle(2); u2->Draw();
  c2->SaveAs(O+"02_elastic_omp_sweep.png");

  // ---------- 3. inelastic: which Ptolemy coupling reproduces FRESCO ----------
  auto*c3=new TCanvas("c3","",1200,520); c3->Divide(2,1);
  const char* on[3]={"INELOCA1","INELOCA2","INELOCA3"};
  int oc[3]={kRed+1,kGreen+2,kMagenta+1};
  auto*gF=rd(F+"p14C_inel_161_7012_L2_dsdo_ex2.dat");
  double fn=gF->Eval(85);
  auto*gFn=new TGraph(); for(int i=0;i<gF->GetN();i++){ double x=gF->GetX()[i];
    if(x>=20&&x<=140) gFn->SetPoint(gFn->GetN(),x,gF->GetY()[i]/fn); }
  c3->cd(1); gPad->SetGridx(); gPad->SetGridy();
  gFn->SetTitle("^{14}C 7.012 MeV 2^{+}: shape, normalised at 85#circ;#theta_{cm} [deg];shape (arb.)");
  gFn->SetLineColor(kBlue+1); gFn->SetLineWidth(4); gFn->SetMinimum(0.5); gFn->SetMaximum(2.6);
  gFn->Draw("AL");
  auto*l3=new TLegend(0.40,0.64,0.92,0.89); l3->SetBorderSize(0); l3->SetFillStyle(0);
  l3->SetTextSize(0.038); l3->AddEntry(gFn,"FRESCO (deformed potential)","l");
  const char* lb[3]={"Ptolemy INELOCA1  (rms 0.087)","Ptolemy INELOCA2  (rms 0.131)","Ptolemy INELOCA3  (rms 0.210)"};
  for(int i=0;i<3;i++){ auto*g=rd(P+TString::Format("inel_7012_%s.dat",on[i]));
    double n=g->Eval(85); auto*q=new TGraph();
    for(int j=0;j<g->GetN();j++){ double x=g->GetX()[j];
      if(x>=20&&x<=140) q->SetPoint(q->GetN(),x,g->GetY()[j]/n); }
    q->SetLineColor(oc[i]); q->SetLineWidth(2); q->SetLineStyle(i==0?1:2); q->Draw("L same");
    l3->AddEntry(q,lb[i],"l"); }
  l3->Draw();
  c3->cd(2); gPad->SetGridx(); gPad->SetGridy();
  for(int i=0;i<3;i++){ auto*g=rd(P+TString::Format("inel_7012_%s.dat",on[i]));
    double n=g->Eval(85); auto*q=new TGraph();
    for(int j=0;j<g->GetN();j++){ double x=g->GetX()[j]; if(x<20||x>140) continue;
      double f=gF->Eval(x)/fn; if(f>0) q->SetPoint(q->GetN(),x,(g->GetY()[j]/n)/f); }
    q->SetLineColor(oc[i]); q->SetLineWidth(2); q->SetLineStyle(i==0?1:2);
    if(i==0){ q->SetTitle("ratio Ptolemy / FRESCO;#theta_{cm} [deg];shape ratio");
              q->SetMinimum(0.5); q->SetMaximum(1.6); q->Draw("AL"); } else q->Draw("L same"); }
  auto*u3=new TLine(20,1,140,1); u3->SetLineColor(kBlack); u3->SetLineStyle(2); u3->Draw();
  TLatex t3; t3.SetNDC(); t3.SetTextSize(0.040);
  t3.DrawLatex(0.15,0.84,"INELOCA1 matches to 2% out to 105#circ");
  t3.DrawLatex(0.15,0.78,"backward tail is a real model difference");
  c3->SaveAs(O+"03_inelastic_coupling_vs_fresco.png");
  printf("\n  wrote 3 figures to %s\n",O.Data());
}
