#include "/home/yassid/Data/C15d_statefit/excm.h"
// Truth check: take the MC deuteron, apply the SAME two-body missing-mass formula the analysis
// uses, and verify it returns the Ex that was generated.  If this fails, nothing downstream means
// anything.
void truth1(const char* f,double exTrue,double Ebeam=191.25){
  TFile fm(f); TTree* t=(TTree*)fm.Get("cbmsim");
  TClonesArray* mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
  Pars P; P.ebeam=Ebeam; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=false;
  const double md=2.0141017781*931.49401;
  TH1D h("h","",200,exTrue-3,exTrue+3); double sth=0,ske=0; long n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    for(int j=0;j<mc->GetEntriesFast();j++){
      auto* m=(FairMCTrack*)mc->At(j);
      if(m->GetPdgCode()!=1000010020) continue;        // the deuteron
      if(m->GetMotherId()!=-1) continue;               // primary only
      double px=m->GetPx()*1000,py=m->GetPy()*1000,pz=m->GetPz()*1000;
      double p=std::sqrt(px*px+py*py+pz*pz), ke=std::sqrt(p*p+md*md)-md;
      double th=std::acos(pz/p)*180.0/M_PI;
      double e=exOf(P,Ebeam,th,ke);
      if(std::isfinite(e)){ h.Fill(e); sth+=th; ske+=ke; n++; }
      break; } }
  printf("  %-28s  N %4ld  <theta_lab> %5.1f  <KE> %5.2f   Ex reco %+.4f +- %.4f  (true %+.2f)\n",
         gSystem->BaseName(f),n,n?sth/n:0,n?ske/n:0,h.GetMean(),h.GetRMS(),exTrue);
}
void truth(){ truth1("/home/yassid/Data/C15dd_sim/tc_0.0_sim.root",0.0);
                 truth1("/home/yassid/Data/C15dd_sim/tc_6.68_sim.root",6.68); }
