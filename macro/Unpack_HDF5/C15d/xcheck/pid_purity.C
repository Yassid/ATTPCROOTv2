// PURITY: of all tracks in a Brho slice, what fraction sits inside a band (within +-2 sigma of a
// fitted p/d/t peak) rather than in the background between and around them?  This is the one
// "separation" measure that survives normalising the statistics -- if ATTPCROOT simply carries
// more junk per track, that is a real difference and not a density artefact.
double betheRel(double pm,double m,double z){
  double E=sqrt(pm*pm+m*m),b=pm/E,g=E/m,b2=b*b; if(b2<1e-6)b2=1e-6;
  const double I=19.2e-6,me=0.511; double a=2*me*b2*g*g/I; if(a<1.)a=1.;
  return z*z/b2*(log(a)-b2); }
struct F{double mu,sg,amp;bool ok;};
F fitAt(TH1D*h,double g,double hw){ F r{0,0,0,false};
  double lo=g*(1-hw),hi=g*(1+hw); if(lo<0.3)lo=0.3;
  if(h->Integral(h->FindBin(lo),h->FindBin(hi))<120) return r;
  TF1 f("f","gaus(0)+pol1(3)",lo,hi);
  f.SetParameters(h->GetBinContent(h->FindBin(g)),g,0.09*g,1,0);
  f.SetParLimits(0,0,1e9);f.SetParLimits(1,lo,hi);f.SetParLimits(2,0.02*g,0.25*g);
  if(h->Fit(&f,"RQ0N")!=0) return r;
  r.mu=f.GetParameter(1);r.sg=fabs(f.GetParameter(2));r.amp=f.GetParameter(0);
  r.ok=r.sg>0&&r.amp>0; return r; }
void purity(){
  const double M[3]={938.272,1875.613,2808.921};
  struct P{const char*f;const char*n;double il,ih,sc;};
  P ps[2]={{"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root","ATTPCROOT",1045,1225,1.0},
           {"/home/yassid/Data/spyral/C15d_a2091/points_spyral.root","Spyral",1000,1250,1.0/1.726}};
  const double blo[5]={0.40,0.50,0.60,0.70,0.85}, bhi[5]={0.50,0.60,0.70,0.85,1.05};
  printf("%-10s %-11s %9s %10s %10s %10s\n","plane","Brho","total","in bands","in-band %","gauss-area %");
  double tp[2]={0,0}, tt[2]={0,0}, ga[2]={0,0};
  for(int ip=0;ip<2;ip++){ auto&p=ps[ip];
    TFile f(p.f); TTree*t=(TTree*)f.Get("pts"); if(!t) continue;
    float x,y,ic;int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    for(int k=0;k<5;k++){
      TH1D*h=new TH1D(Form("P_%s_%d",p.n,k),"",400,0,45); h->SetDirectory(nullptr);
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
      double br=.5*(blo[k]+bhi[k]),pm=299.792458*br,rp=betheRel(pm,M[0],1);
      double tot=h->Integral(), inb=0, garea=0; double muP=pk[0];
      for(int s=0;s<3;s++){
        double pred = s? muP*sqrt(betheRel(pm,M[s],1)/rp) : muP;
        double best=-1,bd=1e9; for(double q:pk){double d=fabs(q-pred)/pred; if(d<bd&&d<0.28){bd=d;best=q;}}
        if(best<0) continue;
        F B=fitAt(h,best,0.17); if(!B.ok) continue;
        inb  += h->Integral(h->FindBin(B.mu-2*B.sg), h->FindBin(B.mu+2*B.sg));
        garea+= B.amp*B.sg*sqrt(2*TMath::Pi())/h->GetBinWidth(1);   // background-subtracted
      }
      printf("%-10s %4.2f-%-6.2f %9.0f %10.0f %9.1f%% %11.1f%%\n",p.n,blo[k],bhi[k],tot,inb,100*inb/tot,100*garea/tot);
      tp[ip]+=inb; tt[ip]+=tot; ga[ip]+=garea;
      delete sm; delete h;
    } }
  printf("\n%-10s %-24s %-24s\n","","overall in-band fraction","overall gaussian-area fraction");
  for(int i=0;i<2;i++) printf("%-10s %-24s %-24s\n",ps[i].n,
     Form("%.1f %%",100*tp[i]/tt[i]), Form("%.1f %%",100*ga[i]/tt[i]));
}
