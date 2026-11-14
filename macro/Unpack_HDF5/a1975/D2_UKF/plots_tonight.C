void plots_tonight(){
  gStyle->SetOptStat(0);
  // 1. Spyral vs genfit, same method, vz 200-400; Spyral also scaled by 1/0.433 to test shape
  std::ifstream in("spyral_same_vz200-400.csv"); std::string l; std::getline(in,l);
  auto gg=new TGraphErrors, gs=new TGraphErrors, gx=new TGraphErrors; double k=291./126;
  while(std::getline(in,l)){ double a,b,yg,ys,dg,eg,ds,es; sscanf(l.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&a,&b,&yg,&ys,&dg,&eg,&ds,&es);
    if(a<15||b>40) continue; double m=0.5*(a+b); int n=gg->GetN();
    gg->SetPoint(n,m-0.4,dg); gg->SetPointError(n,0,eg); gs->SetPoint(n,m,ds); gs->SetPointError(n,0,es); gx->SetPoint(n,m+0.4,ds*k); gx->SetPointError(n,0,es*k); }
  TCanvas c1("c1","",900,650); c1.SetGrid(); c1.SetLeftMargin(0.12);
  auto fr=c1.DrawFrame(13,0,42,8.5); fr->SetTitle("^{16}C(d,p) bound, v_{z} 200-400, identical method;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
  gg->SetMarkerStyle(20); gg->SetMarkerColor(kBlue+1); gg->SetLineColor(kBlue+1); gg->Draw("P same");
  gs->SetMarkerStyle(21); gs->SetMarkerColor(kRed+1); gs->SetLineColor(kRed+1); gs->Draw("P same");
  gx->SetMarkerStyle(25); gx->SetMarkerColor(kRed+1); gx->SetLineColor(kRed+1); gx->SetLineStyle(2); gx->Draw("P same");
  TLegend L(0.14,0.72,0.60,0.88); L.AddEntry(gg,"genfit (ATTPCROOT)","pe"); L.AddEntry(gs,"Spyral-1.0 InterpSolver","pe"); L.AddEntry(gx,"Spyral #times 2.31 (shape test)","pe"); L.Draw();
  c1.SaveAs("spyral_vs_genfit_bound.png");
  // 2. vz profile per theta_cm band: data vs sim reco vs sim generated (from vzfine.C output, re-run into histos)
}
