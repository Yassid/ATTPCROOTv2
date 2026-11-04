// Is "Spyral looks better" a STATISTICS effect?  ATTPCROOT has 1.9x more tracks, so its background
// is 1.9x denser and the gaps between bands fill in, even at identical band width.  Render both
// with the SAME number of entries and the same colour scale.
void equal(){
  gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  // first count how many pass each plane's own IC cut
  long navail[2]={0,0};
  for(int i=0;i<2;i++){ TFile f(ps[i].f); TTree*t=(TTree*)f.Get("pts");
    float x,y,ic; int np; t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);
      if(ic>=ps[i].il&&ic<=ps[i].ih&&np==1&&x>0&&y>0&&x*ps[i].sc<45&&y<1.6) navail[i]++; } }
  long target=std::min(navail[0],navail[1]);
  printf("available after IC cut: ATTPCROOT %ld, Spyral %ld -> using %ld each\n",navail[0],navail[1],target);
  auto* c=new TCanvas("c","",1500,700); c->Divide(2,1);
  TRandom3 rng(12345);
  for(int i=0;i<2;i++){
    TFile f(ps[i].f); TTree*t=(TTree*)f.Get("pts");
    auto* h=new TH2F(Form("e%d",i), Form("%s  (%ld tracks, equal statistics)",ps[i].n,target),300,0,45,300,0,1.6);
    h->SetDirectory(nullptr);
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    double keep=(double)target/navail[i];
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);
      if(!(ic>=ps[i].il&&ic<=ps[i].ih&&np==1&&x>0&&y>0)) continue;
      double xs=x*ps[i].sc; if(xs>=45||y>=1.6) continue;
      if(rng.Uniform()>keep) continue;
      h->Fill(xs,y); }
    h->GetXaxis()->SetTitle("sqrt(dE/dx) [matched scale]"); h->GetYaxis()->SetTitle("B#rho [T m]");
    h->SetMinimum(0.5); h->SetMaximum(200);              // identical colour scale
    c->cd(i+1); gPad->SetLogz(); gPad->SetRightMargin(0.13); h->Draw("colz");
    printf("  %-10s drawn %.0f entries\n", ps[i].n, h->Integral());
  }
  c->SaveAs("/home/yassid/C15d_plane_equalstats.png");
  printf("wrote /home/yassid/C15d_plane_equalstats.png\n");
}
