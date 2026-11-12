#include "binacc.h"
#include "binbuild.C"
// The excess fraction spikes in a narrow vz window in BOTH codes (ATTPCROOT 92-183 mm, Spyral
// 100-167).  What does the Ex spectrum look like THERE compared with the bulk?  That gives the
// background's shape, model-free.
void vzbkg(){
  struct R{const char* nm; double vlo,vhi,blo,bhi;};
  R rr[2]={{"ATTPCROOT",92,183,275,917},{"Spyral",100,167,200,400}};
  for(int ic=0;ic<2;ic++){
    Cfg&k=CFG[ic]; Pars A=parsOf(k); R&r=rr[ic];
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    auto* hS=new TH1D(Form("s%s",k.nm),"",140,-2,12); hS->SetDirectory(nullptr); // spike window
    auto* hB=new TH1D(Form("b%s",k.nm),"",140,-2,12); hB->SetDirectory(nullptr); // bulk
    long ns=0,nb=0;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      if(tcm<10||tcm>20) continue;
      if(vz>=r.vlo&&vz<r.vhi){ hS->Fill(e); ns++; }
      else if(vz>=r.blo&&vz<r.bhi){ hB->Fill(e); nb++; } }
    f->Close();
    // normalise the bulk to the spike on the g.s. peak, then subtract -> the extra component
    double gS=hS->Integral(hS->FindBin(-0.5),hS->FindBin(0.5));
    double gB=hB->Integral(hB->FindBin(-0.5),hB->FindBin(0.5));
    double sc=gB>0?gS/gB:1;
    printf("\n== %s   vz spike [%.0f,%.0f] N=%ld   vs bulk [%.0f,%.0f] N=%ld  (bulk scaled x%.4f on the g.s.)\n",
           k.nm,r.vlo,r.vhi,ns,r.blo,r.bhi,nb,sc);
    printf("   Ex    :"); for(double X=0.5;X<=9.01;X+=0.5) printf("%7.1f",X); printf("\n");
    printf("   spike :"); for(double X=0.5;X<=9.01;X+=0.5) printf("%7.0f",hS->GetBinContent(hS->FindBin(X))); printf("\n");
    printf("   bulk*s:"); for(double X=0.5;X<=9.01;X+=0.5) printf("%7.1f",sc*hB->GetBinContent(hB->FindBin(X))); printf("\n");
    printf("   EXTRA :"); for(double X=0.5;X<=9.01;X+=0.5)
      printf("%7.1f",hS->GetBinContent(hS->FindBin(X))-sc*hB->GetBinContent(hB->FindBin(X))); printf("\n");
    double tot=0; for(int b=1;b<=hS->GetNbinsX();b++){ double x=hS->GetBinCenter(b);
      if(x>1.5&&x<6) tot+=hS->GetBinContent(b)-sc*hB->GetBinContent(b); }
    printf("   integrated extra over Ex 1.5-6 : %.0f events  (%.1f %% of the spike window)\n",tot,100.0*tot/ns);
  }
}
