#include "excm.h"
// Build Ex spectra in theta_CM bins (projectile convention, as in the 15C(p,p') analysis),
// for both codes.  Writes dd_cmbins.root.
struct Cfg{const char* nm; const char* file; double eb,chi2,keLo,keHi,vzHi,icLo,icHi,kcD,kcP,fw;};
Cfg CFG[2]={
 {"ATTPCROOT","/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",
  191.25,0.2,0.5,18,1100,1074,1201,3913.0,-1.5,1.012},
 {"Spyral","/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root",
  195.0,5.0,0.0,20,400,-1,99999,2535.2,14.5,0.548}};
Pars parsOf(const Cfg&k){
  Pars A; A.ebeam=k.eb; A.chi2=k.chi2; A.thLo=0; A.thHi=180; A.keLo=k.keLo; A.keHi=k.keHi;
  A.vzLo=0; A.vzHi=k.vzHi; A.icLo=k.icLo; A.icHi=k.icHi;
  A.kcDenom=k.kcD; A.kcPivot=k.kcP; A.keOff=0; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true; return A; }
// theta_CM binning, matching the (p,p') analysis
const double TCM0=10, TCMSTEP=10, TCM1=80; const int NCM=(int)((TCM1-TCM0)/TCMSTEP);
void binbuild(){
  TFile o("dd_cmbins.root","RECREATE");
  for(auto&k:CFG){
    Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); if(!f||f->IsZombie()){printf("%s MISSING\n",k.nm);continue;}
    TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    std::vector<TH1D*> h(NCM); 
    for(int i=0;i<NCM;i++){ h[i]=new TH1D(Form("%s_cm%d",k.nm,i),
      Form("%s  #theta_{CM} [%.0f,%.0f);E_{x} (MeV);counts / 100 keV",k.nm,TCM0+i*TCMSTEP,TCM0+(i+1)*TCMSTEP),
      280,-5,23); h[i]->SetDirectory(nullptr); }
    auto* hall=new TH1D(Form("%s_all",k.nm),"",280,-5,23); hall->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      hall->Fill(e);
      int b=(int)((tcm-TCM0)/TCMSTEP); if(b>=0&&b<NCM) h[b]->Fill(e); }
    f->Close();
    o.cd(); hall->Write(); for(int i=0;i<NCM;i++) h[i]->Write();
    printf("\n== %s\n   bin  theta_CM      N     Ex of max   counts at 4.6 / 6.9\n",k.nm);
    for(int i=0;i<NCM;i++){
      if(h[i]->GetEntries()<300){ printf("   %2d   %2.0f-%2.0f   %7.0f   (too few)\n",
        i,TCM0+i*TCMSTEP,TCM0+(i+1)*TCMSTEP,h[i]->GetEntries()); continue; }
      int bm=0; double mx=0;
      for(int b=1;b<=h[i]->GetNbinsX();b++){ double x=h[i]->GetBinCenter(b);
        if(x<2.5) continue; if(h[i]->GetBinContent(b)>mx){mx=h[i]->GetBinContent(b);bm=b;} }
      printf("   %2d   %2.0f-%2.0f   %7.0f     %6.2f       %5.0f / %5.0f\n",
        i,TCM0+i*TCMSTEP,TCM0+(i+1)*TCMSTEP,h[i]->GetEntries(),h[i]->GetBinCenter(bm),
        h[i]->GetBinContent(h[i]->FindBin(4.6)),h[i]->GetBinContent(h[i]->FindBin(6.9))); }
  }
  o.Close(); printf("\nwrote dd_cmbins.root\n");
}
