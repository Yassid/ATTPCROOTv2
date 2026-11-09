// The elastic peak is ~2000 sigma and swamps any threshold-based search, so restrict to the
// INELASTIC region and search there on its own scale.
void ppfind2(){
  gStyle->SetOptStat(0);
  TFile f("/home/yassid/Data/C15d_statefit/pp_ex.root");
  TH1D* h0=(TH1D*)f.Get("hpp"); h0->SetDirectory(nullptr);
  auto* h=new TH1D("hin","15C(p,p') inelastic region;E_{x} [MeV];counts / 100 keV",150,1.0,16.0);
  h->SetDirectory(nullptr);
  for(int b=1;b<=h0->GetNbinsX();b++){ double x=h0->GetBinCenter(b);
    if(x>1.0&&x<16.0) h->SetBinContent(h->FindBin(x), h0->GetBinContent(b)); }
  TSpectrum sp(30);
  TH1* bg=sp.Background(h,10,"nodraw");
  int nf=sp.Search(h,1.3,"nodraw",0.12);
  std::vector<double> pk(sp.GetPositionX(),sp.GetPositionX()+nf);
  std::sort(pk.begin(),pk.end());
  printf("%8s %9s %9s %9s %8s\n","Ex","counts","bkg","excess","signif");
  std::vector<double> keep;
  for(double x:pk){ int b=h->FindBin(x);
    double s=h->Integral(b-2,b+2), bk=bg->Integral(b-2,b+2);
    double ex=s-bk, sg=bk>0?ex/std::sqrt(bk):0;
    printf("%8.2f %9.0f %9.0f %9.0f %8.1f%s\n",x,s,bk,ex,sg,sg>5?"  <<":"");
    if(sg>5) keep.push_back(x); }
  printf("\nINELASTIC STATES (>5 sigma):");
  for(double x:keep) printf("  %.2f",x);
  printf("\n");
  auto* c=new TCanvas("c","",1400,700); c->Divide(1,2);
  c->cd(1); h->Draw("hist"); bg->SetLineColor(kRed); bg->Draw("same");
  auto* hs=(TH1D*)h->Clone("hs"); hs->SetDirectory(nullptr); hs->Add(bg,-1);
  c->cd(2); hs->SetTitle("background-subtracted"); hs->Draw("hist");
  c->SaveAs("/home/yassid/C15p_pp_inelastic.png");
  printf("wrote /home/yassid/C15p_pp_inelastic.png\n");
}
