double betheRel(double brho,double A,double z){
  double m=A*931.494, pm=299.792458*brho*z;
  double E=sqrt(pm*pm+m*m),b=pm/E,g=E/m,b2=b*b; if(b2<1e-9)b2=1e-9;
  const double I=19.2e-6,me=0.511; double a=2*me*b2*g*g/I; if(a<1.)a=1.;
  return z*z/b2*(log(a)-b2); }
void leak3(){
  // the saved gate was drawn on the PEAK plane; on the integral plane the scale is 1.86x larger,
  // so it cannot simply be reused -- quantify how far off it now is.
  std::vector<double> gx,gy;
  std::ifstream in("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/deuteron_C15d.json");
  std::string all((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
  size_t p=all.find("\"vertices\"");
  while(true){ size_t a=all.find('[',p+1); if(a==std::string::npos)break;
    size_t b=all.find(']',a); if(b==std::string::npos)break;
    double u,v; if(sscanf(all.substr(a+1,b-a-1).c_str()," %lf , %lf",&u,&v)==2){gx.push_back(u);gy.push_back(v);} p=b; }
  TCutG cut("d",gx.size()); for(size_t i=0;i<gx.size();++i) cut.SetPoint(i,gx[i],gy[i]);
  const char* fs[2]={"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root",
                     "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root"};
  const char* nm[2]={"peak plane (gate drawn here)","INTEGRAL plane (same gate)"};
  for(int k=0;k<2;k++){
    TFile f(fs[k]); TTree* t=(TTree*)f.Get("pts"); if(!t) continue;
    float x,y,ic;int np;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    // anchor the rails on this plane's own proton
    double K0=0;
    { std::vector<double> pr;
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        if(ic<1045||ic>1225||np!=1||x<=0||y<0.55||y>0.65) continue; pr.push_back(x); }
      std::sort(pr.begin(),pr.end());
      auto* hh=new TH1D("hh","",300,0,pr.back()); hh->SetDirectory(nullptr);
      for(double q:pr) hh->Fill(q); hh->Smooth(3);
      int bm=0; double vm=-1;
      for(int b=hh->FindBin(pr.front()); b<=hh->FindBin(pr[pr.size()/2]); ++b)
        if(hh->GetBinContent(b)>vm){vm=hh->GetBinContent(b);bm=b;}
      K0=hh->GetBinCenter(bm)/sqrt(betheRel(0.60,1,1)); delete hh; }
    long n=0,d=0,tt=0,a=0,none=0; const double W=0.136;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(ic<1045||ic>1225||np!=1||x<=0||y<=0) continue;
      if(!cut.IsInside(x,y)) continue; n++;
      bool dD=fabs(x-K0*sqrt(betheRel(y,2,1)))/(K0*sqrt(betheRel(y,2,1)))<W;
      bool dT=fabs(x-K0*sqrt(betheRel(y,3,1)))/(K0*sqrt(betheRel(y,3,1)))<W;
      bool dA=fabs(x-K0*sqrt(betheRel(y,4,2)))/(K0*sqrt(betheRel(y,4,2)))<W;
      if(dD)d++; else if(dT)tt++; else if(dA)a++; else none++; }
    printf("  %-30s K0 %.4f | inside %7ld : d %5.1f%%  t %5.1f%%  alpha %4.1f%%  bkg %5.1f%%\n",
           nm[k], K0, n, 100.*d/n, 100.*tt/n, 100.*a/n, 100.*none/n);
  }
}
