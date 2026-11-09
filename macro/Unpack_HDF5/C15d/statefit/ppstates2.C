#include "exlib.h"
void ppstates2(double eb=198.0){
  gStyle->SetOptStat(0);
  Pars p;
  p.mBeam=15.010599; p.mTarg=1.007825; p.mEject=1.007825; p.mResid=15.010599;
  p.ebeam=eb; p.chi2=5; p.icLo=1037; p.icHi=1252; p.vzLo=-100; p.vzHi=1100;
  TFile f("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/pp/plots_pgate/pp_kin_C15p.root");
  TTree* t=(TTree*)f.Get("dd");
  double ke,th,vz,c2; float ic; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
  auto* h=new TH1D("hpp",Form("15C(p,p')  Ebeam %.1f;E_{x} [MeV];counts / 100 keV",eb),250,-3,22);
  h->SetDirectory(nullptr);
  long n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    if(c2>p.chi2||ic<p.icLo||ic>p.icHi||np!=1) continue;
    if(vz<p.vzLo||vz>p.vzHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double e=exOf(p,p.ebeam,thf,ke);
    if(std::isfinite(e)){ h->Fill(e); n++; } }
  // g.s. position tells us whether Ebeam is right
  TF1 g("g","gaus(0)+pol1(3)",-1.6,1.6);
  g.SetParameters(h->GetMaximum(),0,0.5,10,0); g.SetParLimits(2,0.05,2.0);
  h->Fit(&g,"RQ0N");
  printf("Ebeam %6.1f : %7ld tracks | g.s. mu %+.3f  FWHM %.3f\n",
         eb, n, g.GetParameter(1), 2.355*fabs(g.GetParameter(2)));
  if(fabs(eb-198.0)<0.01){
    auto* c=new TCanvas("c","",1300,650); gPad->SetLogy(); h->Draw("hist");
    c->SaveAs("/home/yassid/C15p_pp_states.png");
    // save for the state hunt
    TFile o("/home/yassid/Data/C15d_statefit/pp_ex.root","RECREATE"); h->Write("hpp"); o.Close();
  }
}
