// The reconstruction efficiency of the C15d chain, as a map in (Ex, theta_CM).
void accfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetPaintTextFormat("4.3f");
  TFile f("/home/yassid/Data/C15dd_sim/acceptance_C15dd.root");
  std::vector<TString> tags; std::vector<double> ex;
  TIter n(f.GetListOfKeys()); TKey* k;
  while((k=(TKey*)n())){ TString s(k->GetName());
    if(!s.BeginsWith("eff_")) continue;
    tags.push_back(s); ex.push_back(TString(s(6,s.Length()-6)).Atof()); }
  std::vector<size_t> idx(tags.size()); for(size_t i=0;i<idx.size();i++) idx[i]=i;
  std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return ex[a]<ex[b];});
  int nx=idx.size();
  auto* h=new TH2D("eff","C15d reconstruction efficiency (PRA level, gain 1.5e5)"
                   ";E_{x} (MeV);#theta_{CM} (deg)",nx,0,nx,6,10,70);
  for(int i=0;i<nx;i++){ TH1D* e=(TH1D*)f.Get(tags[idx[i]]);
    h->GetXaxis()->SetBinLabel(i+1,Form("%.2f",ex[idx[i]]));
    for(int b=1;b<=6;b++) h->SetBinContent(i+1,b,e->GetBinContent(b)); }
  TCanvas* c=new TCanvas("ca","",1150,560);
  gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.12); gPad->SetBottomMargin(0.13);
  h->SetMinimum(0.80); h->SetMaximum(1.0);
  h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleSize(0.05);  h->GetYaxis()->SetTitleSize(0.05);
  h->SetMarkerSize(1.4);
  h->Draw("COLZ TEXT");
  TLatex t; t.SetNDC(true); t.SetTextSize(0.040);
  t.DrawLatex(0.11,0.93,"flat to within 7 % over the whole plane: the efficiency does NOT shape the E_{x} spectrum");
  c->SaveAs("/home/yassid/Data/C15dd_sim/efficiency_C15dd.png");
  printf("wrote efficiency_C15dd.png\n");
}
