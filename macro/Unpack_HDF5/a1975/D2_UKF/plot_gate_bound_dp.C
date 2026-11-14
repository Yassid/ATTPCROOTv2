/// @file plot_gate_bound_dp.C -- bound-peak candidates (Spyral kinematics) on the proton gate, per theta_cm bin.
/// Background: every pattern track (pid/pid_plane_dt_dv1104.root, hpid). Green = inside gate, red = outside.
void plot_gate_bound_dp(const char*fn="effchain_bound_dp_200_400.root",const char*png="gate_bound_dp_200_400.png"){
  gStyle->SetOptStat(0); gStyle->SetPalette(kGreyScale); TColor::InvertPalette();
  TFile f(fn); auto*t=(TNtuple*)f.Get("pd"); auto*G=(TGraph*)f.Get("gate");
  TFile fp("pid/pid_plane_dt_dv1104.root"); auto*hp=(TH2*)fp.Get("hpid");
  const double TB[8]={10,15,20,25,30,35,40,45};
  // frame on the candidates, but keep the gate's low-brho corner visible
  t->Draw("sq:brho","","goff"); double xmx=0,ymx=0; for(int i=0;i<t->GetSelectedRows();++i){ xmx=std::max(xmx,t->GetV1()[i]); ymx=std::max(ymx,t->GetV2()[i]); }
  double X1=std::min(40.,1.15*xmx), Y1=std::min(1.1,1.3*ymx);
  TCanvas c("c","",1700,900); c.Divide(4,2);
  float b,sq,br,in,vz; t->SetBranchAddress("b",&b); t->SetBranchAddress("sq",&sq); t->SetBranchAddress("brho",&br); t->SetBranchAddress("in",&in); t->SetBranchAddress("vz",&vz);
  for(int p=0;p<8;++p){ c.cd(p+1); gPad->SetLogz(); gPad->SetLeftMargin(0.13);
    TH2*h=(TH2*)hp->Clone(Form("h%d",p)); h->GetXaxis()->SetRangeUser(0,X1); h->GetYaxis()->SetRangeUser(0,Y1);
    h->SetTitle(p<7?Form("#theta_{cm} %.0f-%.0f, v_{z} 200-400;#sqrt{dE/dx};B#rho [T m]",TB[p],TB[p+1]):"all bins, v_{z} 200-400;#sqrt{dE/dx};B#rho [T m]"); h->Draw("col");
    G->SetLineColor(kBlue); G->SetLineWidth(2); G->Draw("L same");
    auto gi=new TGraph, go=new TGraph; for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(p<7&&(int)b!=p) continue; (in>0.5?gi:go)->AddPoint(sq,br); }
    gi->SetMarkerStyle(20); gi->SetMarkerSize(0.6); gi->SetMarkerColor(kGreen+2); if(gi->GetN()) gi->Draw("P same");
    go->SetMarkerStyle(20); go->SetMarkerSize(0.8); go->SetMarkerColor(kRed); if(go->GetN()) go->Draw("P same");
    TLatex l; l.SetNDC(); l.SetTextSize(0.05); int n=gi->GetN()+go->GetN();
    l.DrawLatex(0.45,0.84,Form("in %d / %d = %.3f",gi->GetN(),n,n?double(gi->GetN())/n:0)); }
  c.SaveAs(png);
}
