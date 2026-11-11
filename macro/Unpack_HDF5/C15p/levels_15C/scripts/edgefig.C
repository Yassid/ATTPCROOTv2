#include "edgetest.C"
// The kinematic ceiling, and the fact that every theta slice peaks on it.
void edgefig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  double Kb=191.25;
  TCanvas* c=new TCanvas("ce","",1500,650); c->Divide(2,1);
  c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
  auto* g=new TGraph(); for(int i=0;i<=220;i++){ double X=i*0.1; g->SetPoint(i,X,thmax_of(X,Kb)); }
  g->SetTitle("^{15}C(d,d'): the deuteron cannot exceed #theta_{max}(E_{x})"
              ";E_{x} (MeV);#theta_{lab}^{max} (deg)");
  g->SetLineColor(kRed); g->SetLineWidth(3);
  g->GetYaxis()->SetRangeUser(0,95); g->GetXaxis()->SetLimits(0,22);
  g->Draw("AL");
  TBox* b=new TBox(0,70,22,90); b->SetFillColorAlpha(kAzure+1,0.22); b->Draw();
  g->Draw("L same");
  TLatex t; t.SetTextSize(0.034); t.SetTextColor(kAzure+2);
  t.DrawLatex(7.0,80.5,"74 % of the data sits here (#theta 70-90#circ)");
  t.SetTextColor(kRed);
  t.DrawLatex(3.2,58,"#theta_{max} = arcsin(#beta_{d}^{cm}/#beta_{cm})");
  t.SetTextColor(kBlack); t.SetTextSize(0.030);
  t.DrawLatex(3.2,40,Form("#Rightarrow at #theta > 70#circ only E_{x} < %.1f MeV is reachable",ExAtThmax(70.0,Kb)));
  t.DrawLatex(3.2,32,"everything beyond that is NOT ^{15}C(d,d')");
  TLine* l=new TLine(ExAtThmax(70.0,Kb),0,ExAtThmax(70.0,Kb),70);
  l->SetLineStyle(2); l->SetLineColor(kAzure+2); l->SetLineWidth(2); l->Draw();

  c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
  double ts[7]={42,46,50,54,58,62,66}, pk[7]={11.90,10.50,8.70,7.10,5.70,4.50,3.50};
  auto* gm=new TGraph(7); for(int i=0;i<7;i++) gm->SetPoint(i,ts[i],pk[i]);
  gm->SetTitle("where each #theta slice peaks, vs its kinematic edge"
               ";#theta_{lab} (deg);E_{x} of the peak (MeV)");
  gm->SetMarkerStyle(20); gm->SetMarkerSize(1.6); gm->SetMarkerColor(kBlack);
  gm->GetYaxis()->SetRangeUser(0,15); gm->GetXaxis()->SetLimits(38,70);
  gm->Draw("AP");
  auto* ge=new TGraph(); for(int i=0;i<=120;i++){ double th=38+i*(32.0/120); ge->SetPoint(i,th,ExAtThmax(th,Kb)); }
  ge->SetLineColor(kRed); ge->SetLineWidth(3); ge->Draw("L");
  gm->Draw("P same");
  TLegend* lg=new TLegend(0.42,0.74,0.96,0.90); lg->SetTextSize(0.034);
  lg->AddEntry(gm,"observed peak of each 4#circ slice","p");
  lg->AddEntry(ge,"kinematic edge  E_{x}(#theta_{max}=#theta)","l");
  lg->SetBorderSize(0); lg->Draw();
  TLatex t2; t2.SetNDC(true); t2.SetTextSize(0.032);
  t2.DrawLatex(0.17,0.24,"the peak TRACKS the edge in every slice:");
  t2.DrawLatex(0.17,0.18,"a real state would stay at one E_{x}.");
  c->SaveAs("/home/yassid/Data/C15d_statefit/kinematic_edge.png");
  printf("wrote kinematic_edge.png;  Ex reachable at theta>70 is %.2f MeV\n",ExAtThmax(70.0,Kb));
}
