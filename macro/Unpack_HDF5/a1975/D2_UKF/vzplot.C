#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
bool sel(double ke,double th,double vz,double c2,double&ex,double&tc){ if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) return false; return kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc); }
void vzplot(){
  gStyle->SetOptStat(0); double lo[3]={10,20,30},hi[3]={20,30,35}; TH1D*D[3],*S[3],*G[3];
  for(int a=0;a<3;++a){ D[a]=new TH1D(Form("D%d",a),"",20,0,1000); S[a]=new TH1D(Form("S%d",a),"",20,0,1000); G[a]=new TH1D(Form("G%d",a),"",20,0,1000); }
  auto band=[&](double t){ for(int a=0;a<3;++a) if(t>=lo[a]&&t<hi[a]) return a; return -1; };
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int a=band(tc); if(a>=0) D[a]->Fill(vz);} }
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simres_dump.root"); TTree*t=(TTree*)f.Get("r"); double ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("th",&th);t->SetBranchAddress("vz",&vz);t->SetBranchAddress("c2",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int a=band(tc); if(a>=0) S[a]->Fill(vz);} }
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simdz.root"); TTree*t=(TTree*)f.Get("s"); float vz,tc; t->SetBranchAddress("vz",&vz);t->SetBranchAddress("tc",&tc);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); int a=band(tc); if(a>=0) G[a]->Fill(vz);} }
  TCanvas c("c","",1500,500); c.Divide(3,1);
  for(int a=0;a<3;++a){ c.cd(a+1); gPad->SetGrid(); gPad->SetLeftMargin(0.13);
    double m=std::max({D[a]->GetMaximum(),S[a]->GetMaximum(),G[a]->GetMaximum()})*1.15;
    G[a]->SetMaximum(m); G[a]->SetTitle(Form("bound, #theta_{cm} %.0f-%.0f;v_{z} [mm] (reco; generated = true);counts / 50 mm",lo[a],hi[a]));
    G[a]->SetLineColor(kGray+1); G[a]->SetFillColor(kGray); G[a]->SetFillStyle(3004); G[a]->Draw("hist");
    S[a]->SetLineColor(kRed+1); S[a]->SetLineWidth(2); S[a]->Draw("hist same");
    D[a]->SetMarkerStyle(20); D[a]->SetMarkerColor(kBlue+1); D[a]->SetLineColor(kBlue+1); D[a]->Draw("E same");
    TLine l1(200,0,200,m),l2(400,0,400,m),l3(915,0,915,m); for(auto*L:{&l1,&l2}){L->SetLineStyle(2); L->SetLineColor(kGreen+2); L->SetLineWidth(2); L->Draw();} l3.SetLineStyle(3); l3.Draw();
    if(a==0){ auto*L=new TLegend(0.15,0.68,0.62,0.88); L->AddEntry(D[a],"data (genfit)","pe"); L->AddEntry(S[a],"sim reco","l"); L->AddEntry(G[a],"sim generated","f"); L->AddEntry(&l1,"adopted window 200-400","l"); L->AddEntry(&l3,"sim generation ends 915","l"); L->Draw(); } }
  c.SaveAs("vz_profile_bound.png");
}
