#include "resfitC.C"
// THE deliverable that survives: the 6.68 MeV angular distribution.
// Adopted configuration:
//   vz in [300,900]   -- above the Ex calibration break, below the acceptance cliff
//   TIE12 ON          -- the 4.43/4.93 doublet is unresolvable and must be tied
//   resonance set D12 + 5.98 + 6.68 -- 7.38 and 8.28 are NOT determined (yields swing through 0)
//   continuum  p3*(M) q(M) / Ex^n, evaluated at the zero-point-corrected Ex
// Stability is quoted from refitting under a 5 % change in statistics (VZLO 300 -> 350).
void final668(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::ONETERM=true; RH::BKG=RH::HALON; RH::TIE12=true;
  bool keep[6]={true,true,true,true,false,false};
  for(int i=0;i<RH::NR;i++) RH::USE[i]=keep[i];
  Cfg b0=CFG[0], b1=CFG[1];
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",780,580);
  gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetGridy(); gPad->SetLogy();
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  TGraphErrors* g[2]; double ymax=0,ymin=1e30;
  for(int ic=0;ic<2;ic++){
    std::vector<std::vector<double>> A[2];
    for(int pass=0;pass<2;pass++){
      CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO = pass? 350 : 300;
      BB_QUIET=true; binbuild();
      A[pass].clear(); rhRun(CFG[ic],NCM,false,&A[pass]); }
    g[ic]=new TGraphErrors();
    printf("\n==== %s : 6.68 MeV angular distribution (adopted configuration) ====\n",CFG[ic].nm);
    printf("   %-12s %12s %12s   %s\n","theta_CM","yield","refit","stability");
    for(int b=0;b<NCM;b++){
      double y0=(b<(int)A[0].size()&&2<(int)A[0][b].size())?A[0][b][2]:0;
      double y1=(b<(int)A[1].size()&&2<(int)A[1][b].size())?A[1][b][2]:0;
      if(y0<=1e-3) continue;
      double x=TCM0+(b+0.5)*TCMSTEP;
      // uncertainty: the larger of sqrt(N) and the refit spread, scaled out of the 6% stat change
      double sr=std::fabs(y1/0.94-y0);
      double e=std::max(std::sqrt(std::max(1.0,y0)),sr);
      int p=g[ic]->GetN(); g[ic]->SetPoint(p,x,y0); g[ic]->SetPointError(p,TCMSTEP/2,e);
      ymax=std::max(ymax,y0+e); ymin=std::min(ymin,std::max(0.5,y0-e));
      printf("   %4.0f-%-7.0f %12.1f %12.1f   %+6.1f%%\n",
             TCM0+b*TCMSTEP,TCM0+(b+1)*TCMSTEP,y0,y1,100.0*(y1/0.94/y0-1)); }
    g[ic]->SetMarkerStyle(mk[ic]); g[ic]->SetMarkerSize(1.3);
    g[ic]->SetMarkerColor(col[ic]); g[ic]->SetLineColor(col[ic]); g[ic]->SetLineWidth(2);
  }
  auto* fr=gPad->DrawFrame(TCM0-2,0.5*ymin,TCM0+NCM*TCMSTEP+2,2.5*ymax);
  fr->SetTitle("15C(d,d')  6.68 MeV state;#theta_{CM} [deg]  (projectile convention);yield [counts / bin]");
  for(int ic=0;ic<2;ic++) g[ic]->Draw("P");
  auto* lg=new TLegend(0.58,0.74,0.90,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(g[0],"ATTPCROOT","pl"); lg->AddEntry(g[1],"Spyral","pl"); lg->Draw();
  auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.028); tx->SetTextColor(kGray+3);
  tx->DrawLatex(0.17,0.25,"#splitline{v_{z}#in[300,900] mm, doublet tied, continuum p_{3}*q/E_{x}^{n}}"
                          "{errors include the refit spread, not statistics alone}");
  tx->SetTextSize(0.026);
  tx->DrawLatex(0.17,0.185,"the only state whose per-angle yield is reproducible; 7.38 and 8.28 are NOT determined");
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/angdist_668_adopted.png");
  printf("\n   wrote angdist_668_adopted.png\n");
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RH::TIE12=false;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true; BB_QUIET=true; binbuild();
}
