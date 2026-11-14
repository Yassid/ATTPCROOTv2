/// @file countz_bound_dp.C -- bound-peak counts per theta_cm bin vs vertex z (data), with the g.s. sim's
/// expected shape (reco-selected, same cuts) scaled to the data in vz 200-400. Selection = angdist_bound_dp.C.
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
const int NT=7; double TB[NT+1]={10,15,20,25,30,35,40,45};
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
int tb(double t){ for(int i=0;i<NT;++i) if(t>=TB[i]&&t<TB[i+1]) return i; return -1; }
bool sel(double ke,double th,double vz,double c2,double&ex,double&tc){ if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) return false; return kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc); }
void countz_bound_dp(){
  gStyle->SetOptStat(0);
  const int NZ=16; TH1D *hd[NT],*hs[NT]; for(int b=0;b<NT;++b){ hd[b]=new TH1D(Form("d%d",b),Form("#theta_{cm} %.0f-%.0f;v_{z} [mm];bound counts / 50 mm",TB[b],TB[b+1]),NZ,150,950); hs[b]=new TH1D(Form("s%d",b),"",NZ,150,950); }
  TH2D h2("h2","data bound-peak counts;v_{z} [mm];#theta_{cm} [deg]",NZ,150,950,NT,TB);
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int b=tb(tc); if(b<0) continue; hd[b]->Fill(vz); h2.Fill(vz,tc);} }
  { TFile f("/mnt/f/a1975_C16_dp_sim/gs_s4001_simres_dump.root"); TTree*t=(TTree*)f.Get("r"); double ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("th",&th);t->SetBranchAddress("vz",&vz);t->SetBranchAddress("c2",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); double ex,tc; if(!sel(ke,th,vz,c2,ex,tc)||ex<-1.5||ex>0.729) continue; int b=tb(tc); if(b>=0) hs[b]->Fill(vz);} }
  TCanvas c("c","",1700,900); c.Divide(4,2);
  printf("counts per 100 mm (data | sim scaled to 200-400)\n  thcm  "); for(int z=150;z<950;z+=100) printf("  %3d-%3d  ",z,z+100); printf("\n");
  for(int b=0;b<NT;++b){ c.cd(b+1); gPad->SetLeftMargin(0.14);
    double dn=hd[b]->Integral(hd[b]->FindBin(201),hd[b]->FindBin(399)), sn=hs[b]->Integral(hs[b]->FindBin(201),hs[b]->FindBin(399)); if(sn>0) hs[b]->Scale(dn/sn);
    hd[b]->SetMaximum(1.35*std::max(hd[b]->GetMaximum()+sqrt(hd[b]->GetMaximum()),hs[b]->GetMaximum())); hd[b]->SetMinimum(0);
    hd[b]->SetMarkerStyle(20); hd[b]->Draw("E"); hs[b]->SetLineColor(kRed); hs[b]->SetLineWidth(2); hs[b]->Draw("hist same");
    TBox*bx=new TBox(200,0,400,hd[b]->GetMaximum()); bx->SetFillColorAlpha(kGreen+1,0.12); bx->Draw(); hd[b]->Draw("E same");
    TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.17,0.84,Form("N = %.0f",hd[b]->Integral()));
    printf("  %2.0f-%2.0f",TB[b],TB[b+1]); for(int z=0;z<8;++z) printf("  %3.0f | %4.1f",hd[b]->GetBinContent(2*z+1)+hd[b]->GetBinContent(2*z+2),hs[b]->GetBinContent(2*z+1)+hs[b]->GetBinContent(2*z+2)); printf("\n"); }
  c.cd(8); gPad->SetRightMargin(0.14); h2.Draw("colz text");
  c.SaveAs("countz_bound_dp.png");
}
