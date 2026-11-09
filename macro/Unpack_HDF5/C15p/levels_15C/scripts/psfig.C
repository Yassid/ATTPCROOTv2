void psfig(){
  gStyle->SetOptStat(0);
  TFile fd("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TFile ft("/home/yassid/Data/C15d_statefit/ps_alt.root");
  TH1D* hd=(TH1D*)fd.Get("ATTPCROOT"); hd->SetDirectory(nullptr);
  TH1D* hA=(TH1D*)ft.Get("A"); TH1D* hB=(TH1D*)ft.Get("B");
  hA->SetDirectory(nullptr); hB->SetDirectory(nullptr);
  auto* c=new TCanvas("c","",1300,700);
  gPad->SetLogy();
  hd->SetTitle("15C(d,d') continuum vs 3-body phase space (ATTPCROOT);E_{x} [MeV];counts / 100 keV");
  hd->GetXaxis()->SetRangeUser(0,22); hd->SetLineColor(kBlack); hd->SetLineWidth(2);
  hd->GetYaxis()->SetRangeUser(20, hd->GetMaximum()*1.5);
  hd->Draw("hist");
  hA->SetLineColor(kRed);  hA->SetLineWidth(2); hA->SetLineStyle(2); hA->Draw("hist same");
  hB->SetLineColor(kBlue); hB->SetLineWidth(2); hB->Draw("hist same");
  auto* L=new TLegend(0.45,0.68,0.88,0.86); L->SetBorderSize(0);
  L->AddEntry(hd,"data (d,d'), deuteron-gated","l");
  L->AddEntry(hB,"deuteron BREAKUP: 15C + p + n, proton as a deuteron","l");
  L->AddEntry(hA,"15C* -> 14C + n (d + 14C + n)","l");
  L->Draw();
  auto* t=new TLatex(); t->SetTextSize(0.030); t->SetTextColor(kGray+3);
  t->DrawLatex(12.5, hd->GetMaximum()*0.03, "both normalised to the data at E_{x} = 10 MeV");
  c->SaveAs("/home/yassid/C15d_continuum_phasespace.png");
  printf("wrote /home/yassid/C15d_continuum_phasespace.png\n");
}
