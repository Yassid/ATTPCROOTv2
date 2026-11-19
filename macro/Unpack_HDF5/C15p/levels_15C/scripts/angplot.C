TGraphErrors *g00=nullptr,*g01=nullptr,*g10=nullptr,*g11=nullptr;
// angular distributions -- fit from Ex 2.0 (peak excluded from the FIT, background still
// defined from Sn), pure Gaussian peaks, background free per angle, stability-filtered.
void angplot(){ gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",1150,540); c->Divide(2,1);
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  const char* pn[2]={"P1  (4.43 + 4.93)","P2  (5.98 + 6.68)"};
  { double X[]={17.5,22.5,27.5,32.5,42.5,47.5,52.5},Y[]={29.4,20.6,6.7,6.2,2,0.8,1.2},E[]={5.42218,4.53872,3.19362,2.48998,1.41421,0.894427,1.09545};
    g00=new TGraphErrors(7,X,Y,nullptr,E); }
  { double X[]={22.5,27.5,42.5,47.5},Y[]={15.4,9,1,2.6},E[]={3.92428,3,1,1.61245};
    g01=new TGraphErrors(4,X,Y,nullptr,E); }
  { double X[]={12.5,17.5,22.5,32.5,42.5},Y[]={4.7,11.9,22.2,6.4,2.4},E[]={2.16795,3.44964,4.71169,2.52982,1.54919};
    g10=new TGraphErrors(5,X,Y,nullptr,E); }
  { double X[]={17.5,22.5,27.5,37.5,52.5,57.5},Y[]={7.4,19,13,2.7,4,1.1},E[]={2.72029,4.3589,3.60555,1.64317,2,1.04881};
    g11=new TGraphErrors(6,X,Y,nullptr,E); }

  for(int ip=0;ip<2;ip++){
    c->cd(ip+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    gPad->SetGridy(); gPad->SetLogy();
    auto* fr=gPad->DrawFrame(8,0.3,62,120);
    fr->SetTitle(Form("%s;#theta_{CM} [deg]  (projectile convention);yield [counts / 100 keV]",pn[ip]));
    TGraphErrors* gg[2]={ip==0?g00:g10, ip==0?g01:g11};
    for(int ic=0;ic<2;ic++){ if(!gg[ic]) continue;
      gg[ic]->SetMarkerStyle(mk[ic]); gg[ic]->SetMarkerSize(1.5);
      gg[ic]->SetMarkerColor(col[ic]); gg[ic]->SetLineColor(col[ic]); gg[ic]->SetLineWidth(2);
      gg[ic]->Draw("PL"); }
    auto* lg=new TLegend(0.58,0.75,0.92,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
    if(gg[0]) lg->AddEntry(gg[0],"ATTPCROOT","pl");
    if(gg[1]) lg->AddEntry(gg[1],"Spyral","pl");
    lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.027); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.17,0.25,"#splitline{Gaussian peaks; background free per angle, defined}{from S_{n} but FITTED from E_{x}=2 (peak excluded)}");
    tx->DrawLatex(0.17,0.175,"only bins reproducing under a 5% sample change, #chi^{2}_{R}<3");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/angdist_final_gauss.png");
  printf("   wrote angdist_final_gauss.png\n"); }
