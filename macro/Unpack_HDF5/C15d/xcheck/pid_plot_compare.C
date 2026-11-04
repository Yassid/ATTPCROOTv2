// Both planes on MATCHED axes. Spyral's sqrt_dEdx is on a different scale, so it is divided by the
// measured median ratio -- without that the two cannot be compared by eye at all, which is the
// confound in "the Spyral one looks better".
void planeplot(){
  gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
  const char* fa="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root";
  const char* fb="/home/yassid/Data/spyral/C15d_a2091/points_spyral.root";
  // measure the scale ratio from the medians in one Brho slice
  auto med=[&](const char* f,double lo,double hi,double icLo,double icHi){
    TFile ff(f); TTree*t=(TTree*)ff.Get("pts"); float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    std::vector<double> v;
    for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i);
      if(ic<icLo||ic>icHi||np!=1||x<=0||y<lo||y>=hi) continue; v.push_back(x);}
    std::sort(v.begin(),v.end()); return v.empty()?0.0:v[v.size()/2]; };
  double ma=med(fa,0.6,0.8,1045,1225), mb=med(fb,0.6,0.8,1000,1250);
  double k = mb/ma;
  printf("median sqrt_dEdx in Brho 0.6-0.8:  ATTPCROOT %.2f   Spyral %.2f   -> Spyral/ATTPCROOT = %.3f\n", ma, mb, k);
  auto* c=new TCanvas("c","",1500,700); c->Divide(2,1);
  const char* nm[2]={"ATTPCROOT MultiFit  (own scale)","Spyral  (scaled by 1/%.2f to match)"};
  for(int i=0;i<2;i++){
    const char* f = i? fb : fa;
    double icLo = i? 1000:1045, icHi = i? 1250:1225, sc = i? 1.0/k : 1.0;
    TFile ff(f); TTree*t=(TTree*)ff.Get("pts");
    auto* h=new TH2F(Form("h%d",i), i? Form(nm[1],k) : nm[0], 300,0,45, 300,0,1.6);
    h->SetDirectory(nullptr);   // else the open TFile owns it and deletes it on scope exit
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);
      if(ic<icLo||ic>icHi||np!=1) continue; h->Fill(x*sc,y);}
    h->GetXaxis()->SetTitle("sqrt(dE/dx)  [matched scale]");
    h->GetYaxis()->SetTitle("B#rho  [T m]");
    c->cd(i+1); gPad->SetLogz(); gPad->SetRightMargin(0.13); h->Draw("colz");
    printf("  %-34s entries in frame: %.0f\n", h->GetTitle(), h->Integral());
  }
  c->SaveAs("/home/yassid/C15d_plane_compare.png");
  printf("wrote /home/yassid/C15d_plane_compare.png\n");
}
