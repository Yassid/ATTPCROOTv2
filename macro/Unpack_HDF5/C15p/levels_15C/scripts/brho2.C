#include "binacc.h"
#include "binbuild.C"
// a1975 recorded a "bulge" the spiral fit puts in the data deuteron band at Brho 0.22-0.30 T m,
// i.e. slow deuterons.  Brho 0.26 is p = 78 MeV/c = KE 1.6 MeV -- exactly the background
// population here.  Join to the PID points and compare the Brho and track length of the
// background events against the elastic ones.
void brho2(){
  setCM(10,5,60);
  TFile* fp=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root");
  TTree* pt=(TTree*)fp->Get("pts");
  float sq,br,al,dx; int rn,ev,tid,ncl;
  pt->SetBranchAddress("sqrtdedx",&sq); pt->SetBranchAddress("brho",&br);
  pt->SetBranchAddress("arclen",&al);   pt->SetBranchAddress("dedx",&dx);
  pt->SetBranchAddress("run",&rn); pt->SetBranchAddress("event",&ev); pt->SetBranchAddress("trackID",&tid);
  pt->SetBranchAddress("nclusters",&ncl);
  struct P{float br,arc,ncl,lnm;};
  std::unordered_map<long long,P> M; M.reserve(pt->GetEntries()*1.4);
  for(Long64_t i=0;i<pt->GetEntries();i++){ pt->GetEntry(i);
    if(sq<=0||br<=0) continue;
    M[((long long)rn<<40)|((long long)ev<<10)|(tid&0x3FF)]={br,al,(float)ncl,(float)std::log(sq*br)}; }
  fp->Close();
  Cfg& k=CFG[0]; Pars A=parsOf(k);
  TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np,run,event,trackID;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  t->SetBranchAddress("run",&run); t->SetBranchAddress("event",&event); t->SetBranchAddress("trackID",&trackID);
  struct S{const char* n; double lo,hi; long N; double sb,sb2,sa,sn,sl;};
  S w[3]={{"g.s.   Ex -0.5..0.5",-0.5,0.5,0,0,0,0,0,0},
          {"bkg    Ex  3.0..5.0", 3.0,5.0,0,0,0,0,0,0},
          {"bkg    Ex  7.0..9.0", 7.0,9.0,0,0,0,0,0,0}};
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
    double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>20) continue;
    auto it=M.find(((long long)run<<40)|((long long)event<<10)|(trackID&0x3FF));
    if(it==M.end()) continue;
    for(auto&s:w) if(e>=s.lo&&e<s.hi){ s.N++; s.sb+=it->second.br; s.sb2+=it->second.br*it->second.br;
                                        s.sa+=it->second.arc; s.sn+=it->second.ncl; s.sl+=it->second.lnm; } }
  f->Close();
  printf("\n== ATTPCROOT theta_CM bin 0: PID observables of the background vs the elastic\n");
  printf("   %-22s %7s  %16s %9s %8s %10s\n","window","N","<Brho> (T m)","<arclen>","<nclus>","<ln(sqrt(dEdx)Brho)>");
  for(auto&s:w){ if(!s.N) continue;
    double mb=s.sb/s.N, rb=std::sqrt(std::max(0.0,s.sb2/s.N-mb*mb));
    printf("   %-22s %7ld  %7.4f +-%6.4f %9.0f %8.1f %10.4f\n",
           s.n,s.N,mb,rb,s.sa/s.N,s.sn/s.N,s.sl/s.N); }
  printf("\n   a1975 records the spiral-fit 'bulge' in the deuteron band at Brho 0.22-0.30 T m.\n");
}
