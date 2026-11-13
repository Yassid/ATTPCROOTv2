#include "/home/yassid/Data/C15d_statefit/excm.h"
// The +0.21/+0.29 offset is the beam losing energy before it reacts (maxELoss = 10 MeV, so the
// vertex is spread along z and the beam energy there is below the nominal 191.25).  The analysis
// handles this with Ebeam(vz) = ebeam - (ebzSlope/1000)*(vz - ebzPivot).  Apply the SAME
// correction using the MC vertex and the offset should collapse.
void one(const char* f,double exTrue,double slope,double Ebeam=191.25){
  TFile fm(f); TTree* t=(TTree*)fm.Get("cbmsim");
  TClonesArray* mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
  Pars P; P.ebeam=Ebeam; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=true; P.ebzSlope=slope; P.ebzPivot=500;
  const double md=2.0141017781*931.49401;
  TH1D h("h","",400,exTrue-4,exTrue+4); double svz=0; long n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    for(int j=0;j<mc->GetEntriesFast();j++){
      auto* m=(FairMCTrack*)mc->At(j);
      if(m->GetPdgCode()!=1000010020||m->GetMotherId()!=-1) continue;
      double px=m->GetPx()*1000,py=m->GetPy()*1000,pz=m->GetPz()*1000;
      double p=std::sqrt(px*px+py*py+pz*pz), ke=std::sqrt(p*p+md*md)-md;
      double th=std::acos(pz/p)*180.0/M_PI;
      double vz=m->GetStartZ()*10.0;                    // cm -> mm, the analysis unit
      double e=exOf(P,ebAt(P,vz,true),th,ke);
      if(std::isfinite(e)){ h.Fill(e); svz+=vz; n++; }
      break; } }
  printf("   slope %5.1f MeV/m :  <vz> %6.1f mm   Ex reco %+.4f +- %.4f   (true %+.2f, offset %+.4f)\n",
         slope,n?svz/n:0,h.GetMean(),h.GetRMS(),exTrue,h.GetMean()-exTrue);
}
void truth2(){
  // The mean offset is NOT removable by Ebeam(vz): <vz> sits at the pivot, so the correction only
  // takes out the vz-DEPENDENT part.  The offset is the mean beam-energy loss before the reaction
  // (maxELoss = 10 MeV sampled uniformly, so ~5 MeV on average).  The DATA has the same effect and
  // the analysis absorbs it by TUNING Ebeam so the g.s. sits at zero -- 191.25 is an effective
  // value at the mean vertex, not the injection energy.  The simulation must be calibrated the
  // same way; the acceptance and the response are both insensitive to a 0.2 MeV scale offset.
  const char* F[2]={"/home/yassid/Data/C15dd_sim/tc_0.0_sim.root",
                    "/home/yassid/Data/C15dd_sim/tc_6.68_sim.root"};
  double EX[2]={0.0,6.68};
  for(int i=0;i<2;i++){ printf("  Ex true %.2f\n",EX[i]);
    for(double sl : {0.0,11.0,22.0}) one(F[i],EX[i],sl); }
}
