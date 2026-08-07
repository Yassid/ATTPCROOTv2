void accfig(){
  gStyle->SetOptStat(0);
  TFile*a=TFile::Open("/mnt/f/a1954_C14_acc/acceptance_merged_gs.root");
  TFile*b=TFile::Open("/mnt/f/a1954_C14_acc/acceptance_merged_ex1.root");
  TH1D*g=(TH1D*)a->Get("hAcc_gs_sum"),*e=(TH1D*)b->Get("hAcc_ex1_sum");
  TH1D*gg=(TH1D*)a->Get("hGen_gs_sum"),*gr=(TH1D*)a->Get("hRec_gs_sum");
  TH1D*eg=(TH1D*)b->Get("hGen_ex1_sum"),*er=(TH1D*)b->Get("hRec_ex1_sum");
  TCanvas*c=new TCanvas("c","acc",1400,560); c->Divide(2,1);
  c->cd(1);
  g->SetTitle("a1954 ^{14}C(p,p') acceptance, 5 seeds, truth-matched;#theta_{cm} [deg];acceptance");
  g->SetMinimum(0); g->SetMaximum(1.09);
  g->SetMarkerStyle(20); g->SetMarkerColor(kAzure+2); g->SetLineColor(kAzure+2); g->Draw("E1");
  e->SetMarkerStyle(21); e->SetMarkerColor(kOrange+7); e->SetLineColor(kOrange+7); e->Draw("E1 same");
  // shade the recommended usable band
  TBox*ok=new TBox(30,0,110,1.09); ok->SetFillColorAlpha(kGreen-9,0.25); ok->Draw("same");
  TBox*mid=new TBox(110,0,140,1.09); mid->SetFillColorAlpha(kYellow-9,0.25); mid->Draw("same");
  g->Draw("E1 same"); e->Draw("E1 same"); gPad->RedrawAxis();
  TLegend*l=new TLegend(0.40,0.14,0.88,0.34);
  l->AddEntry(g,"g.s.  (E_{x}=0)","lp"); l->AddEntry(e,"1st exc. (E_{x}=6.094 MeV)","lp");
  l->AddEntry(ok,"use directly (#pm2%)","f"); l->AddEntry(mid,"use, #pm7%","f");
  l->SetBorderSize(0); l->Draw();
  c->cd(2);
  gg->SetTitle("generated (truth) vs reconstructed;#theta_{cm} [deg];reactions");
  gg->SetLineColor(kAzure+2); gg->SetLineStyle(2); gg->SetMinimum(0); gg->SetMaximum(1000); gg->Draw("hist");
  gr->SetLineColor(kAzure+2); gr->SetLineWidth(2); gr->Draw("hist same");
  eg->SetLineColor(kOrange+7); eg->SetLineStyle(2); eg->Draw("hist same");
  er->SetLineColor(kOrange+7); er->SetLineWidth(2); er->Draw("hist same");
  TLegend*l2=new TLegend(0.62,0.62,0.93,0.86); l2->SetFillStyle(0);
  l2->AddEntry(gg,"g.s. generated","l"); l2->AddEntry(gr,"g.s. reconstructed","l");
  l2->AddEntry(eg,"ex1 generated","l"); l2->AddEntry(er,"ex1 reconstructed","l");
  l2->SetBorderSize(0); l2->Draw();
  c->SaveAs("/home/yassid/a1954_C14_acceptance_final.png");
  printf("gs  overall %.4f   ex1 overall %.4f\n", gr->Integral()/gg->Integral(), er->Integral()/eg->Integral());
}
