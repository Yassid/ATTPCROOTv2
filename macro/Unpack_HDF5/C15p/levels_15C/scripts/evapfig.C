#include "evap.C"
#include "binbuild.C"
// Why the background is not fusion evaporation: it fails on the angular distribution and on the
// Ex shape, in opposite directions, for every temperature.
void evapfig(){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.05,"t");
  gRandom->SetSeed(11);
  double Kb=191.25, Es=EV::Estar(Kb), T=EV::temperature(Es), B=EV::coulombBarrier();
  double pb=std::sqrt(Kb*Kb+2*Kb*EV::M15C);
  TLorentzVector W(0,0,pb,Kb+EV::M15C+EV::MD); TVector3 bv=W.BoostVector();
  Pars A; A.ebeam=Kb; A.mBeam=15.0105993; A.mTarg=2.014102; A.mEject=2.014102; A.mResid=15.0105993;
  A.kcOn=false; A.ebzOn=false;
  TH1D hth("hth","",9,0,90); TH1D hex("hex","",60,-2,22);
  const int N=8000000; long nl=0;
  for(int i=0;i<N;i++){
    double eps; while(true){ eps=B+gRandom->Exp(T)+gRandom->Exp(T); if(eps<Es) break; }
    double p=std::sqrt(eps*eps+2*eps*EV::MD);
    double c=-1+2*gRandom->Rndm(), st=std::sqrt(std::max(0.0,1-c*c)), ph=2*M_PI*gRandom->Rndm();
    TLorentzVector d(p*st*std::cos(ph),p*st*std::sin(ph),p*c,eps+EV::MD);
    d.Boost(bv);
    double th=d.Theta()*180/M_PI, ke=d.E()-EV::MD;
    hth.Fill(th);
    if(ke<0.5||ke>18||th<14||th>86) continue;
    double e=exOf(A,Kb,th,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,Kb,th,ke);
    if(tcm>=10&&tcm<20){ hex.Fill(e); nl++; } }
  TCanvas* c=new TCanvas("ce","",1600,660); c->Divide(2,1);
  // ---- left: angular ----
  c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14);
  double dat[9]={0.0,1.5,2.8,3.4,4.3,6.5,7.2,35.4,38.8};
  auto* hd=new TH1D("hd",";#theta_{lab} (deg);fraction of the sample (%)",9,0,90);
  auto* he=new TH1D("he","",9,0,90);
  for(int b=1;b<=9;b++){ hd->SetBinContent(b,dat[b-1]);
    he->SetBinContent(b,100.0*hth.GetBinContent(b)/N); }
  hd->SetTitle("angular distribution: evaporation goes FORWARD, the data is BACKWARD"
               ";#theta_{lab} (deg);fraction of the sample (%)");
  hd->SetMaximum(45); hd->SetLineColor(kBlack); hd->SetLineWidth(3); hd->SetFillColorAlpha(kGray,0.5);
  hd->GetXaxis()->SetTitleSize(0.05); hd->GetYaxis()->SetTitleSize(0.05);
  hd->Draw("hist");
  he->SetLineColor(kRed+1); he->SetLineWidth(3); he->Draw("hist same");
  TLegend* l1=new TLegend(0.17,0.74,0.60,0.88); l1->SetBorderSize(0); l1->SetTextSize(0.040);
  l1->AddEntry(hd,"data (deuteron gate)","f"); l1->AddEntry(he,"fusion evaporation","l"); l1->Draw();
  TLatex t; t.SetNDC(true); t.SetTextSize(0.037);
  t.DrawLatex(0.17,0.66,"only 1.3 % of evaporated deuterons");
  t.DrawLatex(0.17,0.61,"reach #theta_{lab} > 70#circ, where 74 % of");
  t.DrawLatex(0.17,0.56,"the data sits.");
  // ---- right: Ex shape ----
  c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14);
  auto* fr=gPad->DrawFrame(2.5,0.1,9.5,4.0);
  fr->SetTitle("E_{x} shape where the background lives (#theta_{CM} 10-20#circ)"
               ";E_{x} (MeV);normalised at E_{x} = 4");
  fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
  gPad->SetLogy();
  double n4=hex.GetBinContent(hex.FindBin(4.0)); if(n4<=0) n4=1;
  auto* ge=new TGraph(); for(double X=3;X<=9.01;X+=0.5)
    ge->SetPoint(ge->GetN(),X,hex.GetBinContent(hex.FindBin(X))/n4);
  ge->SetLineColor(kRed+1); ge->SetLineWidth(4); ge->Draw("L");
  for(int i=0;i<2;i++){ double kap=i?0.24:0.29; auto* g=new TGraph();
    for(double X=3;X<=9.01;X+=0.25) g->SetPoint(g->GetN(),X,std::exp(-kap*(X-4.0)));
    g->SetLineColor(i?kAzure+2:kGreen+2); g->SetLineWidth(4); g->SetLineStyle(i?2:1); g->Draw("L"); }
  TLegend* l2=new TLegend(0.17,0.20,0.72,0.38); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.037);
  TGraph* gd1=new TGraph(); gd1->SetLineColor(kGreen+2); gd1->SetLineWidth(4);
  TGraph* gd2=new TGraph(); gd2->SetLineColor(kAzure+2); gd2->SetLineWidth(4); gd2->SetLineStyle(2);
  l2->AddEntry(ge,"fusion evaporation (T = 4.2 MeV)","l");
  l2->AddEntry(gd1,"fitted background, ATTPCROOT (#kappa = 0.29)","l");
  l2->AddEntry(gd2,"fitted background, Spyral (#kappa = 0.24)","l");
  l2->Draw();
  TLatex t2; t2.SetNDC(true); t2.SetTextSize(0.037);
  t2.DrawLatex(0.17,0.84,"evaporation RISES, the background FALLS.");
  t2.DrawLatex(0.17,0.79,"No temperature from 1 to 20 MeV fixes either:");
  t2.DrawLatex(0.17,0.74,"the slope stays #geq -0.03 against -0.24 .. -0.29.");
  c->SaveAs("/home/yassid/Data/C15d_statefit/evaporation_excluded.png");
  printf("wrote evaporation_excluded.png\n");
}
