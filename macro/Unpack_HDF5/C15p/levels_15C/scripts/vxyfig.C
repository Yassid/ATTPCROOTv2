#include "binacc.h"
#include "binbuild.C"
// Where do the background events START?  Elastic deuterons leave the beam axis; a deuteron from a
// SECOND scatter, tens of cm down the first deuteron's path, does not.  Nor would a mis-measured
// elastic event move its own vertex.  Both codes, theta_CM 10-20 deg.
void vxyfig(){
  setCM(10,5,60);
  TChain ch("kin"); ch.Add("/home/yassid/C15d_fit/*_kin_d.root");
  Int_t run,event,track; Double_t vx,vy;
  ch.SetBranchAddress("run",&run); ch.SetBranchAddress("event",&event);
  ch.SetBranchAddress("trackID",&track);
  ch.SetBranchAddress("vx",&vx); ch.SetBranchAddress("vy",&vy);
  std::unordered_map<long long,std::pair<float,float>> M; M.reserve(ch.GetEntries()*1.4);
  for(Long64_t i=0;i<ch.GetEntries();i++){ ch.GetEntry(i);
    M[((long long)run<<40)|((long long)event<<10)|(track&0x3FF)]={(float)vx,(float)vy}; }
  const char* wn[3]={"elastic peak  |E_{x}| < 0.5","background  E_{x} 3-5","background  E_{x} 7-9"};
  double wl[3]={-0.5,3.0,7.0}, wh[3]={0.5,5.0,9.0};
  int col[3]={kBlack,kRed+1,kAzure+2};
  TCanvas* c=new TCanvas("c","c",1150,520); c->Divide(2,1);
  TLegend* lg=nullptr;
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic]; Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2; float icf; int r2,e2,t2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);t->SetBranchAddress("ic",&icf);
    t->SetBranchAddress("run",&r2);t->SetBranchAddress("event",&e2);t->SetBranchAddress("trackID",&t2);
    std::vector<float> gx,gy; std::vector<std::pair<double,long long>> keep;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icf<A.icLo||icf>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke), eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>20) continue;
      long long key=((long long)r2<<40)|((long long)e2<<10)|(t2&0x3FF);
      auto it=M.find(key); if(it==M.end()) continue;
      if(e>-0.5&&e<0.5){ gx.push_back(it->second.first); gy.push_back(it->second.second); }
      keep.push_back({e,key}); }
    f->Close();
    auto med=[](std::vector<float>&v){ std::nth_element(v.begin(),v.begin()+v.size()/2,v.end()); return v[v.size()/2]; };
    float cx=med(gx), cy=med(gy);
    c->cd(ic+1); gPad->SetLeftMargin(0.13); gPad->SetGridy();
    bool first=true;
    for(int w=0;w<3;w++){
      auto* h=new TH1D(Form("h%d%d",ic,w),"",100,0,100);
      for(auto&p:keep){ if(p.first<wl[w]||p.first>=wh[w]) continue;
        auto&v=M[p.second]; h->Fill(std::hypot(v.first-cx,v.second-cy)); }
      if(h->Integral()<20) continue;
      // cumulative fraction
      auto* g=new TGraph(); double tot=h->Integral(), run=0;
      for(int i=1;i<=h->GetNbinsX();i++){ run+=h->GetBinContent(i);
        g->SetPoint(g->GetN(),h->GetBinLowEdge(i+1),100.0*run/tot); }
      g->SetLineColor(col[w]); g->SetLineWidth(3);
      if(first){ g->SetTitle(Form("%s   #theta_{CM} 10-20#circ;vertex distance from the beam axis [mm];events within r  [%%]",k.nm));
        g->GetXaxis()->SetLimits(0,100); g->SetMinimum(0); g->SetMaximum(105);
        g->Draw("AL"); first=false; } else g->Draw("L");
      if(ic==0){ if(!lg){ lg=new TLegend(0.42,0.14,0.92,0.34); lg->SetBorderSize(0); lg->SetFillStyle(0); }
        lg->AddEntry(g,wn[w],"l"); }
    }
    if(ic==0&&lg) lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.028); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.42,0.57,"the background vertices sit TIGHTER on the");
    tx->DrawLatex(0.42,0.52,"beam axis than the elastic peak itself:");
    tx->DrawLatex(0.42,0.47,"not a second scatter, and not a");
    tx->DrawLatex(0.42,0.42,"mis-measured elastic event either");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/background_vertex_radius.png");
  printf("   wrote background_vertex_radius.png\n");
}
