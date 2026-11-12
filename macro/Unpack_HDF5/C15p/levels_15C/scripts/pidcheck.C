#include "binacc.h"
#include "binbuild.C"
// Are the low-Ex excess events actually DEUTERONS?
// Join dd_kin (kinematics) to points_mfi (PID) on (run, event, trackID) and compare the PID
// variable ln(sqrt(dE/dx) * Brho) -- which is proportional to ln(mass), charge cancelling --
// for the excess window against the g.s. window, in the same theta_CM bin and at the same KE.
void pidcheck(){
  // --- load the PID points into a hash
  TFile* fp=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root");
  TTree* pt=(TTree*)fp->Get("pts");
  float sq,br,al,dx; int rn,ev,tid,ncl;   // these branches are Float_t, not Double_t
  pt->SetBranchAddress("sqrtdedx",&sq); pt->SetBranchAddress("brho",&br);
  pt->SetBranchAddress("arclen",&al);   pt->SetBranchAddress("dedx",&dx);
  pt->SetBranchAddress("run",&rn); pt->SetBranchAddress("event",&ev); pt->SetBranchAddress("trackID",&tid);
  pt->SetBranchAddress("nclusters",&ncl);
  struct P{float lnm,arc,ncl;};
  std::unordered_map<long long,P> M; M.reserve(pt->GetEntries()*1.4);
  for(Long64_t i=0;i<pt->GetEntries();i++){ pt->GetEntry(i);
    if(sq<=0||br<=0) continue;
    long long key=((long long)rn<<40)|((long long)ev<<10)|(tid&0x3FF);
    M[key]={(float)std::log(sq*br),(float)al,(float)ncl}; }
  fp->Close();
  printf("loaded %zu PID points\n",M.size());

  Cfg&k=CFG[0]; Pars A=parsOf(k);
  TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np,run,event,trackID;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  t->SetBranchAddress("run",&run); t->SetBranchAddress("event",&event); t->SetBranchAddress("trackID",&trackID);
  struct S{const char* n; double lo,hi; long N,matched; double s,s2,sa,sc;};
  S w[3]={{"g.s.   Ex -0.5..0.5",-0.5,0.5,0,0,0,0,0,0},
          {"EXCESS Ex  2.5..4.0", 2.5,4.0,0,0,0,0,0,0},
          {"states Ex  6.0..7.5", 6.0,7.5,0,0,0,0,0,0}};
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
    double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>20) continue;
    long long key=((long long)run<<40)|((long long)event<<10)|(trackID&0x3FF);
    auto it=M.find(key);
    for(auto&s:w) if(e>=s.lo&&e<s.hi){ s.N++;
      if(it!=M.end()){ s.matched++; s.s+=it->second.lnm; s.s2+=it->second.lnm*it->second.lnm;
                       s.sa+=it->second.arc; s.sc+=it->second.ncl; } } }
  f->Close();
  printf("\n== ATTPCROOT theta_CM bin 0 : PID of the excess vs the g.s.\n");
  printf("   %-22s %7s %8s  %18s %10s %9s\n","window","N","matched","<ln(sqrt(dEdx)*Brho)>","<arclen>","<nclus>");
  for(auto&s:w){ if(!s.matched) continue;
    double m=s.s/s.matched, r=std::sqrt(std::max(0.0,s.s2/s.matched-m*m));
    printf("   %-22s %7ld %8ld  %9.4f +- %.4f %10.0f %9.1f\n",
           s.n,s.N,s.matched,m,r,s.sa/s.matched,s.sc/s.matched); }
  printf("\n   ln(mass) separation d-p = ln(2) = %.4f ; d-t = %.4f\n",std::log(2.0),std::log(3.0/2.0));
}
