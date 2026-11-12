#include "binacc.h"
#include "binbuild.C"
// Hypothesis: the Ex 2.5-4 excess is the GROUND STATE leaking out, through a residual error in
// the Ebeam(vz) correction.  The g.s. holds 45127 events in bin 0 against 1572 in the excess, so
// a 3.5 % mis-placed tail would do it.
// Test: is the g.s. ridge FLAT in vz after the correction, and does the excess cluster in vz?
void vztest(){
  for(auto&k:CFG){
    Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    auto* h2=new TH2D(Form("v%s",k.nm),"",12,0,A.vzHi,140,-2,12); h2->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      if(tcm<10||tcm>20) continue;
      h2->Fill(vz,e); }
    f->Close();
    printf("\n== %s  theta_CM bin 0:  g.s. centroid and the excess fraction, per vz slice\n",k.nm);
    printf("   vz slice     N      g.s. mu    g.s. FWHM   N(Ex 2.5-4)  / N(g.s.)\n");
    for(int b=1;b<=h2->GetNbinsX();b++){
      TH1D* py=h2->ProjectionY(Form("py%s%d",k.nm,b),b,b);
      if(py->GetEntries()<300) continue;
      TF1 g("g","gaus(0)+pol1(3)",-1.5,1.5);
      g.SetParameters(py->GetMaximum(),0,0.4,2,0); g.SetParLimits(2,0.05,1.5);
      py->Fit(&g,"RQ0N");
      double ngs=py->Integral(py->FindBin(-0.5),py->FindBin(0.5));
      double nex=py->Integral(py->FindBin(2.5),py->FindBin(4.0));
      printf("   %4.0f-%4.0f %7.0f   %+7.3f     %6.3f     %7.0f    %6.2f %%\n",
        h2->GetXaxis()->GetBinLowEdge(b),h2->GetXaxis()->GetBinUpEdge(b),py->GetEntries(),
        g.GetParameter(1),2.355*std::fabs(g.GetParameter(2)),nex,ngs>0?100.0*nex/ngs:0.0); }
  }
}
