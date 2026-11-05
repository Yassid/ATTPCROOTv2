// The a1975 (d,p) method: project ln(sqrt(dE/dx) * Brho). For Z=1, Brho ~ p and dE/dx ~ (m/p)^2,
// so sqrt(dE/dx)*Brho ~ m -- the variable is MOMENTUM-INDEPENDENT and p,d,t sit ln2 = 0.693 and
// ln1.5 = 0.405 apart whatever the Brho. One histogram, all the statistics, known spacing.
void lnproj(){
  gStyle->SetOptStat(0);
  struct P{const char*f;const char*n;double il,ih;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250}};
  auto* c=new TCanvas("c","",1500,900); c->Divide(1,2);
  TH1D* H[2];
  for(int i=0;i<2;i++){
    TFile f(ps[i].f); TTree* t=(TTree*)f.Get("pts");
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x); t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic); t->SetBranchAddress("npulse",&np);
    // normalise each plane by its own median of sqrt(dE/dx)*Brho so the axes coincide
    std::vector<double> v;
    for(Long64_t j=0;j<t->GetEntries();j++){ t->GetEntry(j);
      if(ic<ps[i].il||ic>ps[i].ih||np!=1||x<=0||y<=0) continue; v.push_back(x*y); }
    std::sort(v.begin(),v.end()); double med=v[v.size()/2];
    H[i]=new TH1D(Form("ln%d",i),"",400,-2.5,2.5); H[i]->SetDirectory(nullptr);
    for(double q:v) H[i]->Fill(log(q/med));
    printf("%-10s  %zu tracks,  median sqrt(dE/dx)*Brho = %.3f\n", ps[i].n, v.size(), med);
    H[i]->SetTitle(Form("%s;ln(#sqrt{dE/dx}#upoint B#rho / median)  #propto ln m;tracks",ps[i].n));
    H[i]->SetLineColor(i? kRed+1 : kBlue+1); H[i]->SetLineWidth(2);
    c->cd(i+1); H[i]->Draw("hist");
    // mark the expected spacings from the dominant peak
    int bm=H[i]->GetMaximumBin(); double x0=H[i]->GetBinCenter(bm);
    for(int k=-2;k<=2;k++){
      if(!k) continue;
      double xp = x0 + (k>0? (k==1?log(1.5):log(1.5)+log(4.0/3.0)) : (k==-1?-log(2.0):-log(2.0)-log(1.5)));
      auto* L=new TLine(xp,0,xp,H[i]->GetMaximum()*0.8); L->SetLineColor(kGreen+2); L->SetLineStyle(2); L->Draw();
    }
    printf("            dominant peak at ln = %+.3f\n", x0);
  }
  c->SaveAs("/home/yassid/C15d_lnproj.png");
  printf("wrote /home/yassid/C15d_lnproj.png\n");
}
