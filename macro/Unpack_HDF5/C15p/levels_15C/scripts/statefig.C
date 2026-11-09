double bkgf(double x,double*p){ return std::exp(p[0]+p[1]*x)+p[2]+p[3]*x; }
double m2(double*x,double*p){ double v=bkgf(x[0],p);
  v+=p[4]*std::exp(-0.5*std::pow((x[0]-p[5])/p[6],2));
  v+=p[7]*std::exp(-0.5*std::pow((x[0]-p[8])/p[9],2)); return v; }
void statefig(){
  gStyle->SetOptStat(0); gStyle->SetOptFit(0);
  TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
  const char* nm[2]={"ATTPCROOT","Spyral"};
  double sg0[2]={0.55,0.45};
  auto* c=new TCanvas("c","",1400,800); c->Divide(1,2);
  for(int k=0;k<2;k++){
    TH1D* h=(TH1D*)f.Get(nm[k]); h->SetDirectory(nullptr);
    TF1* fit=new TF1(Form("fit%d",k),m2,2.3,13.5,10);
    fit->SetParameters(std::log(std::max(1.0,(double)h->GetBinContent(h->FindBin(2.3)))),-0.3,
                       h->GetBinContent(h->FindBin(13.5)),0);
    fit->SetParameter(4,h->GetBinContent(h->FindBin(4.7))*0.3); fit->SetParameter(5,4.7); fit->SetParameter(6,sg0[k]);
    fit->SetParameter(7,h->GetBinContent(h->FindBin(6.9))*0.3); fit->SetParameter(8,6.9); fit->SetParameter(9,sg0[k]);
    fit->SetParLimits(5,3.8,5.8); fit->SetParLimits(8,6.2,8.0);
    fit->SetParLimits(6,0.15,1.8); fit->SetParLimits(9,0.15,1.8);
    fit->SetParLimits(4,0,1e7); fit->SetParLimits(7,0,1e7);
    h->Fit(fit,"RQ0N");
    c->cd(k+1); gPad->SetLogy();
    h->GetXaxis()->SetRangeUser(1.0,16.0);
    h->SetTitle(Form("%s  (d,d')  -- fit with the states seen in 15C(p,p');E_{x} [MeV];counts / 100 keV",nm[k]));
    h->SetLineColor(kBlack); h->SetLineWidth(2); h->Draw("hist");
    fit->SetLineColor(kRed); fit->SetLineWidth(2); fit->DrawCopy("same");
    auto* bg=new TF1(Form("bg%d",k),[](double*x,double*p){return bkgf(x[0],p);},2.3,13.5,4);
    for(int i=0;i<4;i++) bg->SetParameter(i,fit->GetParameter(i));
    bg->SetLineColor(kBlue); bg->SetLineStyle(2); bg->DrawCopy("same");
    auto* t=new TLatex(); t->SetTextSize(0.045);
    t->DrawLatex(9.0,h->GetMaximum()*0.25,
      Form("#color[2]{E_{x} = %.3f #pm %.3f}   #color[2]{E_{x} = %.3f #pm %.3f}",
           fit->GetParameter(5),fit->GetParError(5),fit->GetParameter(8),fit->GetParError(8)));
    for(double e:{4.75,7.05}){ auto* l=new TLine(e,0.5,e,h->GetMaximum()*0.12);
      l->SetLineColor(kGreen+2); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw(); }
    t->SetTextSize(0.035); t->SetTextColor(kGreen+2);
    t->DrawLatex(1.3,h->GetMaximum()*0.3,"dotted green = 15C(p,p') states");
  }
  c->SaveAs("/home/yassid/C15d_dd_statefit.png");
  printf("wrote /home/yassid/C15d_dd_statefit.png\n");
}
