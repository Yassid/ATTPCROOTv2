static double exOf(double Eb,double th,double ke,double m3a,double m4a){
  const double u=931.49401; double m1=15.0105993*u,m2=1.00782503*u,m3=m3a*u,m4=m4a*u;
  double E1=Eb+m1,p1=sqrt(E1*E1-m1*m1),E3=ke+m3,p3=sqrt(E3*E3-m3*m3),E4=E1+m2-E3;
  double px=p3*sin(th),pz=p3*cos(th),m4x2=E4*E4-(px*px+(p1-pz)*(p1-pz));
  return (m4x2>0)? sqrt(m4x2)-m4 : -999;
}
/// @file exfit_C15p.C
/// @brief Ex spectrum of a C15p kinematics cache, with the g.s. and the 14C bound cluster fitted.
/// The G.S. IS THE RESOLUTION -- it is the only isolated level here. The 6-7 MeV structure is the
/// unresolved 6.09/6.59/6.73/7.01 multiplet spanning 0.92 MeV, so its width is an envelope.
void exfit_C15p(TString cache = "pp/plots_hdb/pp_kin_C15p_d.root", double Eb = 198.0,
                double m3 = 2.01410178, double m4 = 14.0032420){
  gStyle->SetOptStat(0); gStyle->SetOptFit(0);
  TFile*F=TFile::Open(cache);
  TTree*t=(TTree*)F->Get("dd");
  double ke,th,vz,c2; t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  auto*h=new TH1D("hex","15C(p,d)14C  E_{x};E_{x} [MeV];counts / 100 keV",300,-5,25);
  for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i);
    double ex=exOf(Eb,th*TMath::DegToRad(),ke,m3,m4);
    if(ex>-900) h->Fill(ex);}
  printf("entries %.0f\n",h->GetEntries());
  // ground state
  auto*fg=new TF1("fg","gaus(0)+pol1(3)",-2.0,2.0);
  fg->SetParameters(h->GetMaximum(),0,0.5,10,0);
  h->Fit(fg,"RQ0");
  printf("G.S.    mean %7.3f  sigma %6.3f  FWHM %6.3f MeV   (area %.0f)\n",
     fg->GetParameter(1),fabs(fg->GetParameter(2)),2.3548*fabs(fg->GetParameter(2)),
     fg->GetParameter(0)*fabs(fg->GetParameter(2))*2.5066/0.1);
  // the bound cluster near 6.5
  auto*fc=new TF1("fc","gaus(0)+pol1(3)",4.0,9.5);
  fc->SetParameters(h->GetMaximum()*0.3,6.5,0.8,10,0);
  h->Fit(fc,"RQ0+");
  printf("CLUSTER mean %7.3f  sigma %6.3f  FWHM %6.3f MeV   (area %.0f)\n",
     fc->GetParameter(1),fabs(fc->GetParameter(2)),2.3548*fabs(fc->GetParameter(2)),
     fc->GetParameter(0)*fabs(fc->GetParameter(2))*2.5066/0.1);
  printf("spacing g.s.->cluster = %.3f MeV   (level scheme: 6.09-7.01)\n",
     fc->GetParameter(1)-fg->GetParameter(1));
  auto*c=new TCanvas("c","",1100,700); h->Draw("hist");
  fg->SetLineColor(kRed); fg->Draw("same"); fc->SetLineColor(kGreen+2); fc->Draw("same");
  for(double e : {0.0,6.09,6.59,6.73,7.01,8.18}){auto*l=new TLine(e,0,e,h->GetMaximum()*0.55);
    l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw();}
  c->SaveAs("/home/yassid/pd_ex_spectrum.png");
}
