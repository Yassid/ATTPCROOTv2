TGraphErrors *g00=nullptr,*g01=nullptr,*g10=nullptr,*g11=nullptr;
void gsplot(){ gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",1150,540); c->Divide(2,1);
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  const char* pn[2]={"ground state  (elastic)","5/2^{+}  E_{x} = 0.740 MeV"};
  { double X[]={12.5,17.5,32.5,37.5,42.5,47.5,57.5},Y[]={2172.9,6265.2,47.5,31.6,43.1,27.7,3.7},E[]={73.0064,224.668,6.89202,6.38723,10.866,5.26308,1.92354};
    g00=new TGraphErrors(7,X,Y,nullptr,E); }
  { double X[]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5},Y[]={5034.4,6126,1577.3,67.8,9.9,30.9,31.2,18.7,12.8},E[]={117.06,130.149,39.7152,8.23408,3.19787,6.00638,8.22128,4.32435,3.57771};
    g01=new TGraphErrors(9,X,Y,nullptr,E); }
  { double X[]={12.5,17.5,37.5,42.5,47.5,57.5},Y[]={69.6,217.8,9.1,11.4,12.1,7},E[]={8.34266,18.0128,3.01662,3.63404,3.47851,2.64575};
    g10=new TGraphErrors(6,X,Y,nullptr,E); }
  { double X[]={12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5},Y[]={67.8,96.4,90.8,18.3,3.2,4.6,9.3,6,4.3},E[]={8.23408,9.81835,9.5289,4.27785,1.78885,2.14476,3.04959,2.44949,2.07364};
    g11=new TGraphErrors(9,X,Y,nullptr,E); }

  for(int ip=0;ip<2;ip++){
    c->cd(ip+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    gPad->SetGridy(); gPad->SetLogy();
    auto* fr=gPad->DrawFrame(8,1.0,62, ip? 500.0 : 20000.0);
    fr->SetTitle(Form("%s;#theta_{CM} [deg]  (projectile convention);yield [counts]",pn[ip]));
    TGraphErrors* gg[2]={ip==0?g00:g10, ip==0?g01:g11};
    for(int ci=0;ci<2;ci++){ if(!gg[ci]) continue;
      gg[ci]->SetMarkerStyle(mk[ci]); gg[ci]->SetMarkerSize(1.4);
      gg[ci]->SetMarkerColor(col[ci]); gg[ci]->SetLineColor(col[ci]); gg[ci]->SetLineWidth(2);
      gg[ci]->Draw("PL"); }
    auto* lg=new TLegend(0.58,0.76,0.92,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
    if(gg[0]) lg->AddEntry(gg[0],"ATTPCROOT","pl");
    if(gg[1]) lg->AddEntry(gg[1],"Spyral","pl");
    lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.027); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.17,0.24,"#splitline{fitted Gaussian integrals, g.s. to 10 MeV with no gap}{background free per angle}");
    tx->DrawLatex(0.17,0.175,"NOT acceptance-corrected; only bins stable under a 5% sample change");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/angdist_gs_and_5half.png");
  printf("   wrote angdist_gs_and_5half.png\n"); }
