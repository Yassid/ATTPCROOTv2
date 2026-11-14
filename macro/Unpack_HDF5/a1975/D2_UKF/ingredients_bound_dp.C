#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
/// @file ingredients_bound_dp.C  -- the INGREDIENTS of angdist_bound_dp.C, drawn (same selection code)
/// (copied from angdist_bound_dp.C)
/// 16C(d,p)17C BOUND-state (g.s.+0.217+0.331 unresolved) angular distribution, ABSOLUTE.
/// Built 2026-09-24. Fiducial vz 200-450 mm chosen because data vz shape departs from the g.s. sim
/// beyond 450 (mid-chamber loss, cause OPEN) and 150-200 sits on the reconstruction-onset edge.
/// L_W = 316.98 mb^-1 (full D2 target, beam count x areal density, reference_a1975_luminosity)
///       x 250/971.73 (geometric; uniform gas).
/// Inputs: data  /mnt/f/a1975/caches/.explorer/exp_dp_find805.root (IC-gated slim of dp_kin_find805)
///         sim   /mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root (truth) + gs_s4001_simres_dump.root
///               (per-track reco dump from sim_res_dp.C logic, truth-closest fitted track)
///   root -l -b -q angdist_bound_dp.C
// SIM CONDITIONS both sides: keOff 0, kc off, pivot Ebeam(vz), theta_lab 89-161, KE>1.2, chi2<13.5, Ex -1.5..0.729
// acceptance eps = (sim reco-selected, reco vz in W, reco thcm bin) / (generated, TRUE vz in W, TRUE thcm bin)
const int NT=7; double TB[NT+1]={10,15,20,25,30,35,40,45};
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
bool sel(double ke,double th,double vz,double c2,double&ex,double&tc){
  if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) return false; return kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc); }
int tb(double t){ for(int i=0;i<NT;++i) if(t>=TB[i]&&t<TB[i+1]) return i; return -1; }
void ingredients_bound_dp(double z0W=200,double z1W=450,double Lfull=316.98,double Ltgt=971.73,const char*tag="vz200-450"){
  gSystem->Load("libAtReconstruction.so");
  // generated truth
  std::vector<std::array<double,3>> gen; // tvz, tth, tke
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root"); TTree*t=(TTree*)f.Get("cbmsim"); TClonesArray*mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
    t->SetBranchStatus("*",0); t->SetBranchStatus("MCTrack*",1);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      for(int j=0;j<mc->GetEntriesFast();++j){ auto*m=(AtMCTrack*)mc->At(j); if(!m||m->GetPdgCode()!=2212) continue;
        double px=m->GetPx(),py=m->GetPy(),pz=m->GetPz(),p=sqrt(px*px+py*py+pz*pz)*1000;
        gen.push_back({m->GetStartZ()*10.0, acos(pz*1000/p)*TMath::RadToDeg(), sqrt(p*p+938.272*938.272)-938.272}); break; } } }
  printf("generated protons %zu\n",gen.size());
  // sim reco
  struct S{double ke,th,vz,c2,tke,tth,tvz;}; std::vector<S> sim;
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simres_dump.root"); TTree*t=(TTree*)f.Get("r"); S s; t->SetBranchAddress("ke",&s.ke);t->SetBranchAddress("th",&s.th);t->SetBranchAddress("vz",&s.vz);t->SetBranchAddress("c2",&s.c2);
    t->SetBranchAddress("tke",&s.tke);t->SetBranchAddress("tth",&s.tth);t->SetBranchAddress("tvz",&s.tvz);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); sim.push_back(s);} }
  // data
  struct D{double ke,th,vz,c2;}; std::vector<D> dat;
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); dat.push_back({ke,th,vz,c2});} }
  double z0=z0W,z1=z1W;
  gStyle->SetOptStat(0);
  std::vector<TH1D*> hD,hS;
  for(int b=0;b<NT;++b){ hD.push_back(new TH1D(Form("d%d",b),Form("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts / 0.2 MeV",TB[b],TB[b+1]),45,-5,4));
                         hS.push_back(new TH1D(Form("s%d",b),"",45,-5,4)); }
  TH1D hY("hY",";#theta_{cm} [deg];counts",NT,TB), hB("hB","",NT,TB), hNr("hNr","",NT,TB), hNg("hNg","",NT,TB);
  for(auto&d:dat){ double ex,tc; if(d.vz<z0||d.vz>=z1) continue; if(!sel(d.ke,d.th,d.vz,d.c2,ex,tc)) continue; int b=tb(tc); if(b<0) continue;
    hD[b]->Fill(ex); if(ex>=-1.5&&ex<=0.729) hY.Fill(tc); else if(ex>=-4.0&&ex<-1.5) hB.Fill(tc,2.229/2.5); }
  for(auto&s:sim){ double ex,tc; if(s.vz<z0||s.vz>=z1) continue; if(!sel(s.ke,s.th,s.vz,s.c2,ex,tc)) continue; int b=tb(tc); if(b<0) continue;
    hS[b]->Fill(ex); if(ex>=-1.5&&ex<=0.729) hNr.Fill(tc); }
  for(auto&g:gen){ if(g[0]<z0||g[0]>=z1) continue; double ex,tc; if(!kine2b(Eb(g[0]),g[1]*TMath::DegToRad(),g[2],ex,tc)) continue; hNg.Fill(tc); }
  // (1) Ex per theta bin
  TCanvas c1("c1","",1600,850); c1.Divide(4,2);
  for(int b=0;b<NT;++b){ c1.cd(b+1); gPad->SetLeftMargin(0.14);
    double nd=hD[b]->Integral(hD[b]->FindBin(-1.5),hD[b]->FindBin(0.72)), ns=hS[b]->Integral(hS[b]->FindBin(-1.5),hS[b]->FindBin(0.72));
    if(ns>0) hS[b]->Scale(nd/ns);
    double ym=1.35*std::max(hD[b]->GetMaximum()+sqrt(hD[b]->GetMaximum()),hS[b]->GetMaximum()); if(ym<3) ym=3;
    hD[b]->SetMaximum(ym); hD[b]->SetMinimum(0); hD[b]->SetMarkerStyle(20); hD[b]->SetMarkerSize(0.7); hD[b]->Draw("E");
    TBox*w=new TBox(-1.5,0,0.729,ym); w->SetFillColorAlpha(kGreen+1,0.15); w->Draw();
    TBox*sb=new TBox(-4.0,0,-1.5,ym); sb->SetFillColorAlpha(kOrange,0.15); sb->Draw();
    hS[b]->SetLineColor(kRed); hS[b]->SetLineWidth(2); hS[b]->Draw("hist same"); hD[b]->Draw("E same");
    TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.17,0.84,Form("Y=%.0f  bg=%.1f",hY.GetBinContent(b+1),hB.GetBinContent(b+1)));
    t.SetTextSize(0.045); t.DrawLatex(0.17,0.78,Form("#varepsilon=%.3f",hNg.GetBinContent(b+1)>0?hNr.GetBinContent(b+1)/hNg.GetBinContent(b+1):0)); }
  c1.cd(8); TLegend L(0.05,0.3,0.95,0.8); L.SetBorderSize(0);
  L.AddEntry(hD[0],"data, v_{z} 200-450","pe"); L.AddEntry(hS[0],"sim g.s. (scaled to window)","l");
  TBox bw; bw.SetFillColorAlpha(kGreen+1,0.3); TBox bs; bs.SetFillColorAlpha(kOrange,0.3);
  L.AddEntry(&bw,"yield window -1.5..0.729","f"); L.AddEntry(&bs,"bg sideband -4..-1.5","f"); L.Draw();
  c1.SaveAs(Form("ingredients_ex_per_theta_%s.png",tag));
  // (2) yield, acceptance, result
  TCanvas c2("c2","",1600,500); c2.Divide(3,1);
  c2.cd(1); hY.SetTitle("raw yield in window (black), sideband bg (orange)"); hY.SetLineWidth(2); hY.SetLineColor(kBlack); hY.SetMinimum(0); hY.SetMaximum(1.3*hY.GetMaximum()); hY.Draw("hist E");
  hB.SetLineColor(kOrange+1); hB.SetFillColorAlpha(kOrange,0.4); hB.Draw("hist same");
  c2.cd(2); TGraphAsymmErrors ge; ge.Divide(&hNr,&hNg,"cl=0.683 b(1,1) mode"); ge.SetTitle("acceptance #varepsilon (sim g.s.);#theta_{cm} [deg];#varepsilon");
  ge.SetMarkerStyle(20); ge.SetMarkerColor(kRed+1); ge.SetLineColor(kRed+1); ge.Draw("AP"); ge.GetYaxis()->SetRangeUser(0,1.05); ge.GetXaxis()->SetLimits(10,45);
  c2.cd(3); gPad->SetLogy(); TGraphErrors gr; double LW=Lfull*(z1-z0)/Ltgt;
  for(int b=0;b<NT;++b){ double e=hNg.GetBinContent(b+1)>0?hNr.GetBinContent(b+1)/hNg.GetBinContent(b+1):0; if(e<=0) continue;
    double dO=2*TMath::Pi()*(cos(TB[b]*TMath::DegToRad())-cos(TB[b+1]*TMath::DegToRad())), y=hY.GetBinContent(b+1), bg=hB.GetBinContent(b+1), net=y-bg;
    if(net<=0) continue; int n=gr.GetN(); gr.SetPoint(n,0.5*(TB[b]+TB[b+1]),net/(e*dO*LW)); gr.SetPointError(n,2.5,sqrt(y+bg)/(e*dO*LW)); }
  gr.SetTitle(Form("d#sigma/d#Omega = (Y-bg)/(#varepsilon #Delta#Omega L_{W}),  L_{W}=%.2f mb^{-1};#theta_{cm} [deg];mb/sr",LW));
  gr.SetMarkerStyle(20); gr.SetMarkerColor(kBlue+1); gr.SetLineColor(kBlue+1); gr.Draw("AP");
  c2.SaveAs(Form("ingredients_yield_acc_xs_%s.png",tag));
}
