/// @file ic_ladder_C15d.C
/// @brief IC spectrum and the beam Z ladder for the a2091 D2 runs.
///
///   root -b -q "pid/ic_ladder_C15d.C()"
///
/// At fixed separator Brho every cocktail component has the same velocity, so energy loss goes as
/// Z^2 and the IC peaks form a Z ladder anchored on the tallest (carbon). A peak that does NOT sit
/// at an integer Z is not another element at that velocity -- on these runs there is one at 1365
/// ADC, ratio 1.203 to carbon, Z 6.58, carrying 213k events. It is unidentified and the
/// proton-run window [979.5, 1278.8] cuts straight through it.
// The a2091 beam is a COCKTAIL. At fixed separator Brho every component has the same velocity, so
// energy loss goes as Z^2 and the IC peaks form a Z ladder. Identify it, and check per-run
// stability -- a peak that moves run to run is gain drift, not a species.
void ic_ladder_C15d(){
  gStyle->SetOptStat(0);
  TChain ch("ic"); ch.Add("/home/yassid/C15d_ic/*_ic.root");
  Float_t icmax; Int_t npulse;
  ch.SetBranchAddress("icmax",&icmax); ch.SetBranchAddress("npulse",&npulse);
  auto*hA=new TH1D("hA","a2091 D2: ion chamber;IC max [ADC];events",300,0,3000);
  auto*hS=new TH1D("hS","",300,0,3000);
  long tot=0,sp=0,base=0;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i);
    if(icmax<0) continue; ++tot; hA->Fill(icmax);
    if(icmax<150) ++base;
    if(npulse==1){ ++sp; hS->Fill(icmax); } }
  printf("\n=== IC SPECTRUM, a2091 D2 runs ===\n");
  printf("  %ld entries, %ld single-pulse (%.1f%%), %ld below 150 ADC (%.1f%%)\n",
         tot,sp,100.0*sp/std::max(1L,tot),base,100.0*base/std::max(1L,tot));
  auto*s=(TH1D*)hS->Clone("s"); s->Smooth(2);
  double mx=s->GetMaximum(); int ib=0; double best=0;
  for(int b=1;b<=s->GetNbinsX();++b) if(s->GetBinCenter(b)>800&&s->GetBinContent(b)>best){best=s->GetBinContent(b);ib=b;}
  double aC=s->GetBinCenter(ib);
  printf("\n  structures above 0.2%% of the peak, Z from the tallest = Z 6:\n");
  printf("  %10s %10s %9s %8s %s\n","peak[ADC]","events","ratio/C","-> Z","");
  for(int b=3;b<=s->GetNbinsX()-2;++b){ double v=s->GetBinContent(b);
    if(v<0.002*mx) continue;
    if(!(v>=s->GetBinContent(b-1)&&v>=s->GetBinContent(b+1)&&v>s->GetBinContent(b-2)&&v>s->GetBinContent(b+2))) continue;
    double c=s->GetBinCenter(b); long integ=0;
    for(int k=1;k<=hS->GetNbinsX();++k){double x=hS->GetBinCenter(k); if(x>=c-45&&x<=c+45) integ+=(long)hS->GetBinContent(k);}
    double z=6*sqrt(c/aC);
    printf("  %10.0f %10ld %9.3f %8.2f %s\n",c,integ,c/aC,z,fabs(z-round(z))<0.12?"<--":""); }
  printf("\n  a2091 PROTON runs for comparison: C 1145, O 2055; window [979.5, 1278.8]\n");
  auto*c1=new TCanvas("cic","ic",1100,700); c1->SetLogy();
  hA->SetLineColor(kGray+2); hA->SetLineWidth(2); hA->Draw("hist");
  hS->SetLineColor(kAzure+2); hS->SetLineWidth(2); hS->Draw("hist same");
  for(double g : {979.508,1278.8}){ auto*l=new TLine(g,1,g,hA->GetMaximum());
    l->SetLineColor(kOrange+7); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw(); }
  auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.032);
  tx->SetTextColor(kGray+2);  tx->DrawLatex(0.55,0.86,"all");
  tx->SetTextColor(kAzure+2); tx->DrawLatex(0.55,0.82,"single pulse");
  tx->SetTextColor(kOrange+7);tx->DrawLatex(0.55,0.78,"proton-run window [980,1279]");
  c1->SaveAs("plots/ic_spectrum_a2091_D2.png");
  printf("  wrote plots/ic_spectrum_a2091_D2.png\n");
}
