/// full vz range (0-915, sim generation coverage) vs vz 200-400, plus Pereira-Lopez 2020 bound sum (17.2 MeV/u,
/// DIFFERENT beam energy -- reference only). Filled = use=1 bins.
TGraphErrors* rdcsv(const char*fn,bool all){ auto g=new TGraphErrors; std::ifstream in(fn); std::string l; std::map<std::string,int> c;
  while(std::getline(in,l)){ if(l.empty()||l[0]=='#') continue; std::vector<std::string> t; std::stringstream ss(l); std::string k; while(std::getline(ss,k,',')) t.push_back(k);
    if(l[0]=='t'){ for(size_t i=0;i<t.size();++i) c[t[i]]=i; continue; } if(!all&&std::stod(t[c["use"]])<0.5) continue; double v=std::stod(t[c["dsdo_mb_sr"]]); if(v<=0) continue;
    int n=g->GetN(); g->SetPoint(n,std::stod(t[c["thcm_mid"]]),v); g->SetPointError(n,0.5*(std::stod(t[c["thcm_hi"]])-std::stod(t[c["thcm_lo"]])),std::stod(t[c["err_mb_sr"]])); } return g; }
void plot_compare_windows_dp(){
  gStyle->SetOptStat(0); TCanvas c("c","",950,680); c.SetLogy(); c.SetGrid(); c.SetLeftMargin(0.12);
  TH1F*fr=c.DrawFrame(0,0.05,50,100); fr->SetTitle("^{16}C(d,p)^{17}C bound states (g.s.+0.217+0.331);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
  auto gp=new TGraphErrors("/home/yassid/fair_install/16C_dp/RCS/dp_AngularDistro_boundstates_Eex_no56_7May2018_CM.txt","%lg %lg %lg");
  gp->SetMarkerStyle(25); gp->SetMarkerColor(kGray+2); gp->SetLineColor(kGray+2); gp->Draw("P");
  auto g2=rdcsv("dp_bound_angdist_vz200-400.csv",true); g2->SetMarkerStyle(24); g2->SetMarkerSize(1.2); g2->SetMarkerColor(kBlue+1); g2->SetLineColor(kBlue+1); g2->Draw("P");
  auto gf=rdcsv("dp_bound_angdist_vz0-915.csv",true); gf->SetMarkerStyle(20); gf->SetMarkerSize(1.3); gf->SetMarkerColor(kRed+1); gf->SetLineColor(kRed+1); gf->SetLineWidth(2); gf->Draw("P");
  TLegend L(0.45,0.72,0.89,0.89); L.SetTextSize(0.03);
  L.AddEntry(gf,"this work, v_{z} 0-915 (full, L = 298.5 mb^{-1})","pe"); L.AddEntry(g2,"this work, v_{z} 200-400 (L = 65.2 mb^{-1})","pe");
  L.AddEntry(gp,"Pereira-Lopez 2020, 17.2 MeV/u (reference)","pe"); L.Draw();
  c.SaveAs("dp_bound_angdist_fullrange_vs_window.png"); }
