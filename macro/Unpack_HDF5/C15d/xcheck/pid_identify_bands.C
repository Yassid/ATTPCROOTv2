// IDENTIFY the bands before measuring them. List every significant peak in a Brho slice and
// compare against the Bethe prediction for p, d, t, 3He and alpha, anchored on the proton.
// sqrt_dEdx scales as sqrt(z^2 * f(beta)), and beta at fixed Brho depends on m/z.
double betheRel(double p_MeVc,double m,double z){
  double E=sqrt(p_MeVc*p_MeVc+m*m), beta=p_MeVc/E, g=E/m, b2=beta*beta; if(b2<1e-6)b2=1e-6;
  const double I=19.2e-6, me=0.511; double arg=2*me*b2*g*g/I; if(arg<1.0)arg=1.0;
  return z*z/b2*(log(arg)-b2);
}
void ident(){
  const char* nm[5]={"p","d","t","3He","a"};
  const double M[5]={938.272,1875.613,2808.921,2808.391,3727.379};
  const double Z[5]={1,1,1,2,2};
  struct P { const char* f; const char* name; double icLo,icHi,sc; };
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  const double blo[2]={0.45,0.60}, bhi[2]={0.55,0.75};
  for(auto&p:ps){
    TFile f(p.f); TTree*t=(TTree*)f.Get("pts"); if(!t) continue;
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(int k=0;k<2;k++){
      TH1D* h=new TH1D(Form("i_%s_%d",p.name,k),"",360,0,45); h->SetDirectory(nullptr);
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        if(ic<p.icLo||ic>p.icHi||np!=1||x<=0||y<blo[k]||y>=bhi[k]) continue; h->Fill(x*p.sc); }
      TH1D* sm=(TH1D*)h->Clone("sm"); sm->SetDirectory(nullptr); sm->Smooth(4);
      std::vector<std::pair<double,double>> pk;
      for(int b=5;b<=sm->GetNbinsX()-5;b++){ double c=sm->GetBinContent(b);
        if(c<sm->GetMaximum()*0.02) continue;
        bool ispk=true; for(int d=-4;d<=4;d++) if(d&&sm->GetBinContent(b+d)>c) ispk=false;
        if(ispk) pk.push_back({sm->GetBinCenter(b),c}); }
      double br=0.5*(blo[k]+bhi[k]), pmom=299.792458*br;
      printf("\n== %s  Brho %.2f-%.2f  (%.0f tracks)\n", p.name, blo[k], bhi[k], h->GetEntries());
      printf("   peaks found at x =");
      for(auto&q:pk) printf(" %.1f(%.0f)", q.first, q.second);
      printf("\n");
      if(!pk.empty()){
        double muP=pk[0].first, rp=betheRel(pmom,M[0],Z[0]);
        printf("   anchored on proton = %.2f -> predictions:", muP);
        for(int s=1;s<5;s++) printf("  %s %.1f", nm[s], muP*sqrt(betheRel(pmom,M[s],Z[s])/rp));
        printf("\n");
      }
      delete sm; delete h;
    }
  }
}
