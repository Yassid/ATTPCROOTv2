// angular distributions, final configuration -- only bins that pass the stability test
TGraphErrors *g00=nullptr,*g01=nullptr,*g10=nullptr,*g11=nullptr;
void angplot(){ gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",1150,540); c->Divide(2,1);
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  const char* pn[2]={"P1  (4.43 + 4.93)","P2  (5.98 + 6.68)"};
  { double X[]={17.5,22.5,27.5,32.5}; double Y[]={44.7,55.5,6.2,3.3}; double E[]={6.68581,7.44983,2.48998,1.81659}; int N=4;
    g00 = N? new TGraphErrors(N,X,Y,nullptr,E) : nullptr; }
  { double X[]={17.5,22.5,27.5,47.5}; double Y[]={33.5,18.4,3,1.4}; double E[]={9.7766,4.28952,1.73205,1.18322}; int N=4;
    g01 = N? new TGraphErrors(N,X,Y,nullptr,E) : nullptr; }
  { double X[]={12.5,17.5,22.5,32.5,42.5}; double Y[]={4.2,14.6,46.1,4.6,1.8}; double E[]={2.04939,3.82099,6.7897,2.14476,1.34164}; int N=5;
    g10 = N? new TGraphErrors(N,X,Y,nullptr,E) : nullptr; }
  { double X[]={12.5,17.5,22.5,27.5,37.5}; double Y[]={2.9,10.1,21.6,8.8,1}; double E[]={1.70294,3.17805,4.64758,2.96648,1}; int N=5;
    g11 = N? new TGraphErrors(N,X,Y,nullptr,E) : nullptr; }

  for(int ip=0;ip<2;ip++){
    c->cd(ip+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    gPad->SetGridy(); gPad->SetLogy();
    auto* fr=gPad->DrawFrame(8,0.4,62,200);
    fr->SetTitle(Form("%s;#theta_{CM} [deg]  (projectile convention);yield [counts / 100 keV]",pn[ip]));
    TGraphErrors* gg[2]={ip==0?g00:g10, ip==0?g01:g11};
    for(int ic=0;ic<2;ic++){ if(!gg[ic]) continue;
      gg[ic]->SetMarkerStyle(mk[ic]); gg[ic]->SetMarkerSize(1.5);
      gg[ic]->SetMarkerColor(col[ic]); gg[ic]->SetLineColor(col[ic]); gg[ic]->SetLineWidth(2);
      gg[ic]->Draw("PL"); }
    auto* lg=new TLegend(0.58,0.74,0.92,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
    if(gg[0]) lg->AddEntry(gg[0],"ATTPCROOT","pl");
    if(gg[1]) lg->AddEntry(gg[1],"Spyral","pl");
    lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.028); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.18,0.23,"#splitline{Gaussian peaks, background free per angle,}{fitted g.s. to 10 MeV with no gap}");
    tx->DrawLatex(0.18,0.165,"only bins reproducing under a 5% sample change");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/angdist_final_gauss.png");
  printf("   wrote angdist_final_gauss.png\n");
}
