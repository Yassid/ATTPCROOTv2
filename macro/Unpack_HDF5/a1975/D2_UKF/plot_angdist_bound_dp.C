/// Plot of dp_bound_angdist_vz200-450.csv (angdist_bound_dp.C). Filled = recommended, open = use=0.
void plot_angdist_bound_dp(const char*csv="dp_bound_angdist_pertheta.csv",const char*png="dp_bound_angdist_pertheta.png",bool logy=true,const char*ttl="per-#theta v_{z} ranges",const char*lw="L per bin"){
  std::ifstream in(csv); std::string l; auto g1=new TGraphErrors, g0=new TGraphErrors;
  std::map<std::string,int> col;  // columns read BY NAME so both csv layouts work
  while(std::getline(in,l)){ if(l.empty()||l[0]=='#') continue;
    std::vector<std::string> tk; std::stringstream ss(l); std::string tok; while(std::getline(ss,tok,',')) tk.push_back(tok);
    if(l[0]=='t'){ for(size_t k=0;k<tk.size();++k) col[tk[k]]=k; continue; }
    auto v=[&](const char*n){ return std::stod(tk[col[n]]); };
    auto g=v("use")>0.5?g1:g0; int n=g->GetN(); g->SetPoint(n,v("thcm_mid"),v("dsdo_mb_sr")); g->SetPointError(n,0.5*(v("thcm_hi")-v("thcm_lo")),v("err_mb_sr")); }
  gStyle->SetOptStat(0); TCanvas c("c","",900,650); c.SetLeftMargin(0.12); c.SetGrid(); c.SetLogy(logy);
  double ymax=0, ymin=1e9; for(auto g:{g1,g0}) for(int i=0;i<g->GetN();++i){ ymax=std::max(ymax,g->GetY()[i]+g->GetEY()[i]); ymin=std::min(ymin,g->GetY()[i]); }
  TH1F*fr=logy? c.DrawFrame(8,0.3*ymin,47,15*ymax) : c.DrawFrame(8,0,47,1.15*ymax);
  fr->SetTitle(Form("^{16}C(d,p)^{17}C bound states (g.s.+0.217+0.331), %s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]",ttl));
  g1->SetMarkerStyle(20); g1->SetMarkerSize(1.3); g1->SetMarkerColor(kBlue+1); g1->SetLineColor(kBlue+1); g1->SetLineWidth(2); g1->Draw("P same");
  g0->SetMarkerStyle(24); g0->SetMarkerSize(1.3); g0->SetMarkerColor(kGray+2); g0->SetLineColor(kGray+2); g0->Draw("P same");
  TLegend L(0.14,0.76,0.62,0.88); L.AddEntry(g1,Form("data, stat. errors (%s)",lw),"pe"); L.AddEntry(g0,"excluded (not adopted)","pe"); L.Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.028); t.DrawLatex(0.14,0.71,"no keOff, no #theta corr., E_{beam}(v_{z}) pivot; no systematics");
  c.SaveAs(png);
}
