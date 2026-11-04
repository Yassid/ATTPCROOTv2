// Zoom on the TRITON region, equal statistics, identical colour scale, with the Bethe-predicted
// p/d/t loci overlaid so the third band can be identified by eye rather than guessed at.
double betheRel(double pm,double m,double z){
  double E=sqrt(pm*pm+m*m),b=pm/E,g=E/m,b2=b*b; if(b2<1e-6)b2=1e-6;
  const double I=19.2e-6,me=0.511; double a=2*me*b2*g*g/I; if(a<1.)a=1.;
  return z*z/b2*(log(a)-b2); }
void triton(){
  gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  long nav[2]={0,0};
  for(int i=0;i<2;i++){ TFile f(ps[i].f); TTree*t=(TTree*)f.Get("pts"); float x,y,ic;int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);
      if(ic>=ps[i].il&&ic<=ps[i].ih&&np==1&&x>0&&y>0&&x*ps[i].sc<35&&y>0.25&&y<1.0) nav[i]++;} }
  long tgt=std::min(nav[0],nav[1]);
  printf("triton-region tracks: ATTPCROOT %ld, Spyral %ld -> %ld each\n",nav[0],nav[1],tgt);
  // proton reference: calibrate the locus normalisation from the ATTPCROOT proton band
  const double M[3]={938.272,1875.613,2808.921};
  auto* c=new TCanvas("c","",1500,700); c->Divide(2,1);
  TRandom3 rng(999);
  for(int i=0;i<2;i++){
    TFile f(ps[i].f); TTree*t=(TTree*)f.Get("pts");
    auto*h=new TH2F(Form("t%d",i),Form("%s  triton region (equal stats)",ps[i].n),260,0,35,220,0.25,1.0);
    h->SetDirectory(nullptr);
    float x,y,ic;int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    double keep=(double)tgt/nav[i];
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);
      if(!(ic>=ps[i].il&&ic<=ps[i].ih&&np==1&&x>0&&y>0)) continue;
      double xs=x*ps[i].sc; if(xs>=35||y<=0.25||y>=1.0) continue;
      if(rng.Uniform()>keep) continue; h->Fill(xs,y); }
    h->GetXaxis()->SetTitle("sqrt(dE/dx) [matched]"); h->GetYaxis()->SetTitle("B#rho [T m]");
    h->SetMinimum(0.5); h->SetMaximum(150);
    c->cd(i+1); gPad->SetLogz(); gPad->SetRightMargin(0.13); h->Draw("colz");
    // overlay the three loci, normalised so the proton curve passes through the measured band
    double K = 6.05 / sqrt(betheRel(299.792458*0.65, M[0], 1));   // ATTPCROOT proton at Brho 0.65
    for(int s=0;s<3;s++){
      auto* g=new TGraph();
      for(double br=0.26; br<1.0; br+=0.01){
        double v = K*sqrt(betheRel(299.792458*br, M[s], 1));
        if(v>0&&v<35) g->SetPoint(g->GetN(), v, br); }
      g->SetLineColor(s==0?kWhite:(s==1?kRed:kMagenta)); g->SetLineWidth(2);
      g->SetLineStyle(s==2?2:1); g->Draw("L same"); }
    auto* l=new TLatex(); l->SetNDC(); l->SetTextSize(0.035);
    l->SetTextColor(kWhite);   l->DrawLatex(0.60,0.85,"white = proton");
    l->SetTextColor(kRed);     l->DrawLatex(0.60,0.80,"red = deuteron");
    l->SetTextColor(kMagenta); l->DrawLatex(0.60,0.75,"magenta = triton");
  }
  c->SaveAs("/home/yassid/C15d_triton_region.png");
  printf("wrote /home/yassid/C15d_triton_region.png\n");
}
