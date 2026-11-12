#include "binacc.h"
#include "binbuild.C"
// Where does the low-Ex excess live in the RAW observables?  If it is instrumental (short,
// low-energy tracks) it sits at a fixed low deuteron KE and therefore at DIFFERENT Ex in
// different theta_CM bins.  If it is nuclear it sits at fixed Ex.
void keorigin(){
  for(auto&k:CFG){
    Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    // KE spectrum of events in the excess window (Ex 2.5-4.0) vs a control window (Ex 5.0-6.0)
    auto* hx=new TH1D(Form("x%s",k.nm),"",60,0,15); hx->SetDirectory(nullptr);
    auto* hc=new TH1D(Form("c%s",k.nm),"",60,0,15); hc->SetDirectory(nullptr);
    // and the Ex of the excess per theta_CM bin
    const int NB2=3; double exlo[NB2]={0},exhi[NB2]={0};
    auto* hke=new TH2D(Form("k%s",k.nm),"",6,0,60,60,0,15); hke->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      if(e>2.5&&e<4.0) hx->Fill(ke);
      if(e>5.0&&e<6.0) hc->Fill(ke);
      hke->Fill(tcm,ke); }
    f->Close();
    printf("\n== %s   deuteron KE of events in the excess window vs a control window\n",k.nm);
    printf("   KE      :"); for(double K=0.5;K<=8.01;K+=0.5) printf("%6.1f",K); printf("\n");
    printf("   Ex 2.5-4:"); for(double K=0.5;K<=8.01;K+=0.5) printf("%6.0f",hx->GetBinContent(hx->FindBin(K))); printf("\n");
    printf("   Ex 5.0-6:"); for(double K=0.5;K<=8.01;K+=0.5) printf("%6.0f",hc->GetBinContent(hc->FindBin(K))); printf("\n");
    printf("   <KE> excess %.2f   <KE> control %.2f   MeV\n",hx->GetMean(),hc->GetMean());
    // what KE does Ex=3.2 correspond to, per theta_CM bin?  (two-body, exact)
    printf("   Ex=3.2 maps to deuteron KE:");
    for(int ib=0;ib<5;ib++){ double cm=TCM0+ib*TCMSTEP+TCMSTEP/2, thl,kk;
      labOf(A,k.eb,3.2,cm,thl,kk);
      printf("  cm%2.0f: %.1f MeV (lab %.0f deg)",cm,kk,thl); }
    printf("\n");
  }
}
