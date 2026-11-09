// Refined: smooth falling continuum (exp + linear) + two Gaussians. Three passes per dataset --
// (a) no peaks, for the null chi2; (b) positions FIXED at the (p,p') values; (c) positions FREE,
// to test whether the (d,d') data actually want the (p,p') energies.
double bkgf(double x,double*p){ return std::exp(p[0]+p[1]*x) + p[2] + p[3]*x; }
double m0(double*x,double*p){ return bkgf(x[0],p); }
double m2(double*x,double*p){
  double v=bkgf(x[0],p);
  v += p[4]*std::exp(-0.5*std::pow((x[0]-p[5])/p[6],2));
  v += p[7]*std::exp(-0.5*std::pow((x[0]-p[8])/p[9],2));
  return v; }
void run(const char* nm,double lo,double hi,double sg0){
  TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
  TH1D* h=(TH1D*)f.Get(nm); h->SetDirectory(nullptr);
  double bw=h->GetBinWidth(1);
  // (a) background only
  TF1 f0("f0",m0,lo,hi,4);
  f0.SetParameters(std::log(std::max(1.0,h->GetBinContent(h->FindBin(lo)))),-0.3,
                   h->GetBinContent(h->FindBin(hi)),0);
  h->Fit(&f0,"RQ0N");
  double chi0=f0.GetChisquare(); int nd0=f0.GetNDF();
  // (b) fixed positions
  TF1 fb("fb",m2,lo,hi,10);
  for(int i=0;i<4;i++) fb.SetParameter(i,f0.GetParameter(i));
  fb.SetParameter(4,h->GetBinContent(h->FindBin(4.75))*0.3); fb.FixParameter(5,4.75); fb.SetParameter(6,sg0);
  fb.SetParameter(7,h->GetBinContent(h->FindBin(7.05))*0.3); fb.FixParameter(8,7.05); fb.SetParameter(9,sg0);
  fb.SetParLimits(6,0.15,1.8); fb.SetParLimits(9,0.15,1.8);
  fb.SetParLimits(4,0,1e7); fb.SetParLimits(7,0,1e7);
  h->Fit(&fb,"RQ0N");
  // (c) free positions
  TF1 fc("fc",m2,lo,hi,10);
  for(int i=0;i<10;i++) fc.SetParameter(i,fb.GetParameter(i));
  fc.ReleaseParameter(5); fc.SetParLimits(5,3.8,5.8);
  fc.ReleaseParameter(8); fc.SetParLimits(8,6.2,8.0);
  fc.SetParLimits(6,0.15,1.8); fc.SetParLimits(9,0.15,1.8);
  fc.SetParLimits(4,0,1e7); fc.SetParLimits(7,0,1e7);
  h->Fit(&fc,"RQ0N");
  printf("\n== %-10s  fit %.1f-%.1f   bkg-only chi2/ndf %.2f  |  +2 peaks fixed %.2f  |  free %.2f\n",
         nm, lo, hi, chi0/std::max(1,nd0),
         fb.GetChisquare()/std::max(1,fb.GetNDF()), fc.GetChisquare()/std::max(1,fc.GetNDF()));
  printf("   delta-chi2 vs background-only: fixed %.0f (4 dof), free %.0f (6 dof)\n",
         chi0-fb.GetChisquare(), chi0-fc.GetChisquare());
  const char* lbl[2]={"state A","state B"};
  double ppv[2]={4.75,7.05};
  for(int i=0;i<2;i++){
    int a=4+3*i;
    double A=fc.GetParameter(a), mu=fc.GetParameter(a+1), s=fabs(fc.GetParameter(a+2));
    double eA=fc.GetParError(a), emu=fc.GetParError(a+1);
    double Y=A*s*std::sqrt(2*M_PI)/bw, eY=eA*s*std::sqrt(2*M_PI)/bw;
    printf("   %-8s (p,p') %.2f -> FREE fit %.3f +- %.3f  FWHM %.3f  yield %.0f +- %.0f  shift %+.3f\n",
           lbl[i], ppv[i], mu, emu, 2.355*s, Y, eY, mu-ppv[i]);
  }
}
void statefit2(){ run("ATTPCROOT",2.3,13.5,0.55); run("Spyral",2.3,13.5,0.45); }
