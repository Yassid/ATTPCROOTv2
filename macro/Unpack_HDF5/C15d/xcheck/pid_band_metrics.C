// FINAL per-band comparison. For each Brho slice: find all peaks, identify p/d/t by matching the
// Bethe-predicted positions (anchored on the proton), then fit each at its OBSERVED peak.
// Reports sigma/mu per band and R_dt = (mu_t-mu_d)/(sig_d+sig_t), the deuteron-triton separation.
double betheRel(double pm,double m,double z){
  double E=sqrt(pm*pm+m*m),b=pm/E,g=E/m,b2=b*b; if(b2<1e-6)b2=1e-6;
  const double I=19.2e-6,me=0.511; double a=2*me*b2*g*g/I; if(a<1.)a=1.;
  return z*z/b2*(log(a)-b2);
}
struct F { double mu,sg,amp,bkg; bool ok; };
F fitAt(TH1D*h,double g,double hw){
  F r{0,0,0,0,false}; double lo=g*(1-hw),hi=g*(1+hw); if(lo<0.3)lo=0.3;
  if(h->Integral(h->FindBin(lo),h->FindBin(hi))<120) return r;
  TF1 f("f","gaus(0)+pol1(3)",lo,hi);
  f.SetParameters(h->GetBinContent(h->FindBin(g)),g,0.09*g,1,0);
  f.SetParLimits(0,0,1e9); f.SetParLimits(1,lo,hi); f.SetParLimits(2,0.02*g,0.25*g);
  if(h->Fit(&f,"RQ0N")!=0) return r;
  r.mu=f.GetParameter(1); r.sg=fabs(f.GetParameter(2)); r.amp=f.GetParameter(0);
  r.bkg=std::max(0.0,f.GetParameter(3)+f.GetParameter(4)*r.mu);
  r.ok=r.sg>0&&r.amp>0&&r.mu>lo*1.02&&r.mu<hi*0.98;
  return r;
}
void final(){
  const double M[3]={938.272,1875.613,2808.921};
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  const double blo[5]={0.40,0.50,0.60,0.70,0.85}, bhi[5]={0.50,0.60,0.70,0.85,1.05};
  printf("%-10s %-11s | %-13s | %-19s | %-19s | %7s\n","plane","Brho","p mu s/m","d  mu sig s/m S/B","t  mu sig s/m S/B","R_dt");
  double sumRd[2]={0,0}; int nRd[2]={0,0}; double sumWd[2]={0,0}; int nWd[2]={0,0};
  for(int ip=0; ip<2; ++ip){ auto&p=ps[ip];
    TFile f(p.f); TTree*t=(TTree*)f.Get("pts"); if(!t) continue;
    float x,y,ic; int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(int k=0;k<5;k++){
      TH1D*h=new TH1D(Form("F_%s_%d",p.n,k),"",400,0,45); h->SetDirectory(nullptr);
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        if(ic<p.il||ic>p.ih||np!=1||x<=0||y<blo[k]||y>=bhi[k]) continue; h->Fill(x*p.sc); }
      if(h->GetEntries()<4000){ delete h; continue; }
      TH1D*sm=(TH1D*)h->Clone("s"); sm->SetDirectory(nullptr); sm->Smooth(4);
      std::vector<double> pk;
      for(int b=5;b<=sm->GetNbinsX()-5;b++){ double c=sm->GetBinContent(b);
        if(c<sm->GetMaximum()*0.015) continue; bool q=true;
        for(int d=-4;d<=4;d++) if(d&&sm->GetBinContent(b+d)>c) q=false;
        if(q) pk.push_back(sm->GetBinCenter(b)); }
      if(pk.empty()){ delete sm; delete h; continue; }
      double br=.5*(blo[k]+bhi[k]), pm=299.792458*br, rp=betheRel(pm,M[0],1);
      double muP=pk[0];
      F B[3]; double pred[3];
      for(int s=0;s<3;s++){
        pred[s]= s? muP*sqrt(betheRel(pm,M[s],1)/rp) : muP;
        double best=-1,bd=1e9;
        for(double q:pk){ double d=fabs(q-pred[s])/pred[s]; if(d<bd&&d<0.28){bd=d;best=q;} }
        B[s]= best>0 ? fitAt(h,best,0.17) : F{0,0,0,0,false};
      }
      printf("%-10s %4.2f-%-6.2f |",p.n,blo[k],bhi[k]);
      if(B[0].ok) printf(" %5.2f %5.3f |",B[0].mu,B[0].sg/B[0].mu); else printf(" %5s %5s |","-","-");
      for(int s=1;s<3;s++){ if(B[s].ok) printf(" %5.2f %5.2f %5.3f %4.1f |",B[s].mu,B[s].sg,B[s].sg/B[s].mu,B[s].bkg>0?B[s].amp/B[s].bkg:99.9);
                            else printf(" %5s %5s %5s %4s |","-","-","-","-"); }
      if(B[1].ok&&B[2].ok&&B[2].mu>B[1].mu){ double R=(B[2].mu-B[1].mu)/(B[1].sg+B[2].sg);
        printf(" %7.2f",R); sumRd[ip]+=R; nRd[ip]++; } else printf(" %7s","-");
      if(B[1].ok){ sumWd[ip]+=B[1].sg/B[1].mu; nWd[ip]++; }
      printf("\n"); delete sm; delete h;
    }
  }
  printf("\n%-10s %-28s %-28s\n","","mean deuteron width sigma/mu","mean d-t separation R_dt");
  for(int i=0;i<2;i++)
    printf("%-10s %-28s %-28s\n", ps[i].n,
      nWd[i]? Form("%.4f  (%d slices)", sumWd[i]/nWd[i], nWd[i]) : "-",
      nRd[i]? Form("%.2f  (%d slices)", sumRd[i]/nRd[i], nRd[i]) : "-");
}
