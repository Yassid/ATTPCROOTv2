#include "resfit.C"
// FINAL angular distributions: Ex 3-9, 5 deg bins to 60 deg, five populated resonances with the
// unresolvable 4.43/4.93 doublet tied.  Only the two structures both codes agree on are drawn.
void finalfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.042,"t");
  setCM(10,5,60);
  RF::USE[0]=RF::USE[1]=RF::USE[2]=RF::USE[3]=RF::USE[4]=true; RF::USE[5]=false; RF::TIE12=true;
  double A[2][4][10]; int nb=10;
  for(int ic=0;ic<2;ic++){
    std::vector<std::vector<double>> amp;
    rfRun(CFG[ic],nb,false,&amp);
    for(int b=0;b<nb;b++){ double s=std::sin((TCM0+(b+0.5)*TCMSTEP)*M_PI/180.0);
      for(int j=0;j<4;j++) A[ic][j][b]=(b<(int)amp.size()&&j<(int)amp[b].size())?amp[b][j]/s:0; } }
  RF::TIE12=false;
  TCanvas* c=new TCanvas("cf","",1020,730);
  gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.08);
  auto* fr=gPad->DrawFrame(8,0.004,62,12.0);
  fr->SetTitle("^{15}C(d,d') angular distributions  --  E_{x} 3-9 MeV, 5#circ bins"
               ";#theta_{CM} (deg);d#sigma/d#Omega   (normalised to each maximum)");
  fr->GetXaxis()->SetTitleSize(0.044); fr->GetYaxis()->SetTitleSize(0.044);
  fr->GetYaxis()->SetTitleOffset(1.3);
  int col[2]={kGreen+2,kOrange+7}; int jj[2]={0,2};
  const char* nm[2]={"4.43 + 4.93 MeV (tied doublet)","6.68 MeV"};
  TLegend* lg=new TLegend(0.40,0.74,0.96,0.91);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
  for(int j=0;j<2;j++)
    for(int ic=0;ic<2;ic++){
      double mx=0; for(int b=0;b<nb;b++) mx=std::max(mx,A[ic][jj[j]][b]);
      if(mx<=0) continue;
      auto* g=new TGraph();
      for(int b=0;b<nb;b++){ double v=A[ic][jj[j]][b]/mx;
        if(v>0.004) g->SetPoint(g->GetN(),TCM0+(b+0.5)*TCMSTEP,v); }
      g->SetLineColor(col[j]); g->SetMarkerColor(col[j]); g->SetLineWidth(3);
      g->SetLineStyle(ic?2:1); g->SetMarkerStyle(ic?25:20); g->SetMarkerSize(1.7);
      g->Draw("LP");
      lg->AddEntry(g,Form("%s  %s",nm[j],ic?"(Spyral)":"(ATTPCROOT)"),"lp"); }
  lg->Draw();
  TLatex t; t.SetNDC(true); t.SetTextSize(0.031);
  t.DrawLatex(0.17,0.26,"#color[418]{4.4+4.9} peaks at 18#circ, #color[807]{6.68} at 22#circ -- in BOTH codes.");
  t.DrawLatex(0.17,0.21,"Different shapes #Rightarrow different angular-momentum transfer.");
  t.SetTextSize(0.027);
  t.DrawLatex(0.17,0.155,"10 angular points. The low-E_{x} region (g.s., 5/2^{+}, continuum turn-on,");
  t.DrawLatex(0.17,0.115,"low-energy background) is OUTSIDE the fit window, so none of it enters here.");
  c->SaveAs("/home/yassid/Data/C15d_statefit/angdist_final.png");
  printf("wrote angdist_final.png\n");
}
