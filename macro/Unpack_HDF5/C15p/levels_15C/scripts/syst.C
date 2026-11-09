double bkgf(double x,double*p){ return std::exp(p[0]+p[1]*x)+p[2]+p[3]*x; }
double m2(double*x,double*p){ double v=bkgf(x[0],p);
  v+=p[4]*std::exp(-0.5*std::pow((x[0]-p[5])/p[6],2));
  v+=p[7]*std::exp(-0.5*std::pow((x[0]-p[8])/p[9],2)); return v; }
void one(const char* nm,double sg0){
  TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* h=(TH1D*)f.Get(nm); h->SetDirectory(nullptr);
  printf("\n== %s : fit-range systematic\n   %-12s %10s %10s %9s %9s\n",nm,"range","Ex(A)","Ex(B)","FWHM A","FWHM B");
  std::vector<double> A,B;
  for(double lo : {2.0,2.3,2.6,3.0}) for(double hi : {12.0,13.5,15.0}){
    TF1 fit("fit",m2,lo,hi,10);
    fit.SetParameters(std::log(std::max(1.0,(double)h->GetBinContent(h->FindBin(lo)))),-0.3,
                      h->GetBinContent(h->FindBin(hi)),0);
    fit.SetParameter(4,h->GetBinContent(h->FindBin(4.7))*0.3); fit.SetParameter(5,4.7); fit.SetParameter(6,sg0);
    fit.SetParameter(7,h->GetBinContent(h->FindBin(6.9))*0.3); fit.SetParameter(8,6.9); fit.SetParameter(9,sg0);
    fit.SetParLimits(5,3.8,5.8); fit.SetParLimits(8,6.2,8.0);
    fit.SetParLimits(6,0.15,1.8); fit.SetParLimits(9,0.15,1.8);
    fit.SetParLimits(4,0,1e7); fit.SetParLimits(7,0,1e7);
    h->Fit(&fit,"RQ0N");
    printf("   %4.1f-%-7.1f %10.3f %10.3f %9.3f %9.3f\n",lo,hi,
           fit.GetParameter(5),fit.GetParameter(8),
           2.355*fabs(fit.GetParameter(6)),2.355*fabs(fit.GetParameter(9)));
    A.push_back(fit.GetParameter(5)); B.push_back(fit.GetParameter(8));
  }
  auto stat=[](std::vector<double>&v,const char* l){
    double m=0; for(double q:v)m+=q; m/=v.size();
    double s=0; for(double q:v)s+=(q-m)*(q-m); s=std::sqrt(s/std::max<size_t>(1,v.size()-1));
    printf("   %s mean %.3f  spread(rms) %.3f\n",l,m,s); };
  stat(A,"state A:"); stat(B,"state B:");
}
void syst(){ one("ATTPCROOT",0.55); one("Spyral",0.45); }
