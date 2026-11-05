// Assumption-free: locate the local maxima and minima of the ln-mass projection in each code and
// report the peak-to-valley contrast between consecutive structures. The idealised ln(m) spacing
// assumes dE/dx ~ 1/beta^2 exactly; the Bethe log term shifts the real positions, so imposing the
// windows can land them off-peak.
void findpk(){
  struct P{const char*f;const char*n;double il,ih;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250}};
  for(auto&p:ps){
    TFile f(p.f); TTree* t=(TTree*)f.Get("pts");
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    std::vector<double> v;
    for(Long64_t j=0;j<t->GetEntries();j++){ t->GetEntry(j);
      if(ic<p.il||ic>p.ih||np!=1||x<=0||y<=0) continue; v.push_back(x*y); }
    std::sort(v.begin(),v.end()); double med=v[v.size()/2];
    auto* h=new TH1D(Form("fp_%s",p.n),"",320,-1.6,1.6); h->SetDirectory(nullptr);
    for(double q:v) h->Fill(log(q/med));
    h->Smooth(2);
    std::vector<int> pk, vl;
    for(int b=4;b<=h->GetNbinsX()-4;b++){
      double c=h->GetBinContent(b);
      if(c<h->GetMaximum()*0.015) continue;
      bool ip=true, iv=true;
      for(int d=-3;d<=3;d++){ if(!d) continue;
        if(h->GetBinContent(b+d)>c) ip=false;
        if(h->GetBinContent(b+d)<c) iv=false; }
      if(ip) pk.push_back(b); else if(iv) vl.push_back(b);
    }
    printf("\n== %-10s %zu tracks\n", p.n, v.size());
    printf("   local maxima at ln =");
    for(int b:pk) printf("  %+.3f(%.0f)", h->GetBinCenter(b), h->GetBinContent(b));
    printf("\n   spacings from the dominant peak:");
    int bd=h->GetMaximumBin();
    for(int b:pk) if(b!=bd) printf("  %+.3f", h->GetBinCenter(b)-h->GetBinCenter(bd));
    printf("\n   reference: ln(1/2)=-0.693  ln(3/2)=+0.405  ln(4/2)=+0.693  ln(3)=+1.099\n");
    // contrast between consecutive maxima
    for(size_t i=0;i+1<pk.size();++i){
      double vmin=1e18; int bmin=pk[i];
      for(int b=pk[i]; b<=pk[i+1]; ++b) if(h->GetBinContent(b)<vmin){vmin=h->GetBinContent(b);bmin=b;}
      double lo=h->GetBinContent(pk[i]), hi=h->GetBinContent(pk[i+1]);
      printf("   %+.3f -> %+.3f : valley %.0f at %+.3f, contrast valley/min(peaks) = %.3f\n",
             h->GetBinCenter(pk[i]), h->GetBinCenter(pk[i+1]), vmin, h->GetBinCenter(bmin),
             vmin/std::min(lo,hi));
    }
    delete h;
  }
}
