// Full-statistics check of the claim made from 13k tracks: does the integral plane separate the
// bands more than the peak plane, and does it match Spyral?  ln(sqrt(dE/dx)*Brho) ~ mass, so the
// spacings from the deuteron are known: p -0.578, t +0.330, 3He +0.455, 4He +0.693 (full Bethe).
void finalcmp(){
  struct P{const char*f;const char*n;double il,ih;};
  P ps[3]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT peak",1045,1225},
           {"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root","ATTPCROOT INTEGRAL",1045,1225},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250}};
  for(auto&p:ps){
    TFile f(p.f); TTree* t=(TTree*)f.Get("pts"); if(!t){printf("%-20s MISSING\n",p.n);continue;}
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    std::vector<double> v;
    for(Long64_t j=0;j<t->GetEntries();j++){ t->GetEntry(j);
      if(ic<p.il||ic>p.ih||np!=1||x<=0||y<=0) continue; v.push_back(x*y); }
    std::sort(v.begin(),v.end()); double med=v[v.size()/2];
    auto* h=new TH1D(Form("F_%s",p.n),"",400,-1.6,1.6); h->SetDirectory(nullptr);
    for(double q:v) h->Fill(log(q/med));
    h->Smooth(2);
    int bd=h->GetMaximumBin();
    std::vector<double> sp;
    for(int b=5;b<=h->GetNbinsX()-5;b++){ double c=h->GetBinContent(b);
      if(c<h->GetMaximum()*0.02) continue; bool ip=true;
      for(int d=-4;d<=4;d++) if(d&&h->GetBinContent(b+d)>c) ip=false;
      if(ip&&b!=bd) sp.push_back(h->GetBinCenter(b)-h->GetBinCenter(bd)); }
    printf("%-20s %8zu tracks | spacings from d:", p.n, v.size());
    for(double d:sp) if(fabs(d)<1.25) printf(" %+.3f", d);
    printf("\n");
    // d and t widths -> resolving power
    TF1 fd("fd","gaus(0)+pol1(3)",-0.26,0.22);
    fd.SetParameters(h->GetMaximum(),0,0.10,10,0); fd.SetParLimits(2,0.02,0.30);
    h->Fit(&fd,"RQ0N");
    double tg=0; for(double d:sp) if(d>0.25&&d<0.70) { tg=d; break; }
    if(tg>0){
      TF1 ft("ft","gaus(0)+pol1(3)",tg-0.15,tg+0.15);
      ft.SetParameters(h->GetBinContent(h->FindBin(tg)),tg,0.07,10,0);
      ft.SetParLimits(1,tg-0.08,tg+0.08); ft.SetParLimits(2,0.02,0.22);
      h->Fit(&ft,"RQ0N");
      double sd=fabs(fd.GetParameter(2)), st=fabs(ft.GetParameter(2));
      printf("%-20s %8s | d sigma %.4f  t mu %+.3f sigma %.4f  ->  R = %.2f\n",
             "", "", sd, ft.GetParameter(1), st, (ft.GetParameter(1)-fd.GetParameter(1))/(sd+st));
    }
    delete h;
  }
  printf("\n  expected from full Bethe:  p -0.578   t +0.330   3He +0.455   4He +0.693\n");
}
