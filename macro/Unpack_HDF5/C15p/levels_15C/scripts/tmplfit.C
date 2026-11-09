// Fit the (d,d') continuum with a PHYSICAL template instead of exp+linear, and check whether the
// extracted state energies depend on that choice (a background-model systematic).
TH1D* gTmpl=nullptr;
double withTmpl(double*x,double*p){
  double b = (gTmpl? gTmpl->Interpolate(x[0]) : 0.0)*p[0];
  b += p[1]*std::exp(-0.5*std::pow((x[0]-p[2])/p[3],2));
  b += p[4]*std::exp(-0.5*std::pow((x[0]-p[5])/p[6],2));
  return b; }
double expLin(double*x,double*p){
  double b=std::exp(p[0]+p[1]*x[0])+p[2]+p[3]*x[0];
  b += p[4]*std::exp(-0.5*std::pow((x[0]-p[5])/p[6],2));
  b += p[7]*std::exp(-0.5*std::pow((x[0]-p[8])/p[9],2));
  return b; }
void run(const char* dataName,const char* tname){
  TFile fd("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* h=(TH1D*)fd.Get(dataName); h->SetDirectory(nullptr);
  TFile ft("/home/yassid/Data/C15d_statefit/ps_alt.root");
  gTmpl=(TH1D*)ft.Get(tname); gTmpl->SetDirectory(nullptr);
  double lo=2.3,hi=20.0;
  TF1 f("f",withTmpl,lo,hi,7);
  f.SetParameters(1.0,h->GetBinContent(h->FindBin(4.7))*0.3,4.7,0.8,
                     h->GetBinContent(h->FindBin(6.9))*0.3,6.9,0.5);
  f.SetParLimits(0,0,1e4); f.SetParLimits(1,0,1e8); f.SetParLimits(4,0,1e8);
  f.SetParLimits(2,3.8,5.8); f.SetParLimits(5,6.2,8.0);
  f.SetParLimits(3,0.15,1.8); f.SetParLimits(6,0.15,1.8);
  h->Fit(&f,"RQ0N");
  printf("   %-10s tmpl %-18s chi2/ndf %6.2f | A %.3f +- %.3f  B %.3f +- %.3f\n",
    dataName,tname,f.GetChisquare()/std::max(1,f.GetNDF()),
    f.GetParameter(2),f.GetParError(2),f.GetParameter(5),f.GetParError(5));
}
void runEmp(const char* dataName){
  TFile fd("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* h=(TH1D*)fd.Get(dataName); h->SetDirectory(nullptr);
  double lo=2.3,hi=20.0;
  TF1 f("f2",expLin,lo,hi,10);
  f.SetParameters(std::log(std::max(1.0,(double)h->GetBinContent(h->FindBin(lo)))),-0.3,
                  h->GetBinContent(h->FindBin(hi)),0,
                  h->GetBinContent(h->FindBin(4.7))*0.3,4.7,0.8,
                  h->GetBinContent(h->FindBin(6.9))*0.3,6.9,0.5);
  f.SetParLimits(5,3.8,5.8); f.SetParLimits(8,6.2,8.0);
  f.SetParLimits(6,0.15,1.8); f.SetParLimits(9,0.15,1.8);
  f.SetParLimits(4,0,1e8); f.SetParLimits(7,0,1e8);
  h->Fit(&f,"RQ0N");
  printf("   %-10s tmpl %-18s chi2/ndf %6.2f | A %.3f +- %.3f  B %.3f +- %.3f\n",
    dataName,"exp+linear (empirical)",f.GetChisquare()/std::max(1,f.GetNDF()),
    f.GetParameter(5),f.GetParError(5),f.GetParameter(8),f.GetParError(8));
}
void tmplfit(){
  printf("background model comparison (fit 2.3-20.0 MeV):\n");
  for(const char* d : {"ATTPCROOT","Spyral"}){
    run(d,"B"); run(d,"A"); runEmp(d); printf("\n"); }
}
