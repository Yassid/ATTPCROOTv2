#include "/home/yassid/Data/C15d_statefit/excm.h"
// The simulation knows the true beam-energy loss, so it can VALIDATE the analysis' Ebeam(vz)
// slope from first principles: the correct slope is the one that minimises the Ex spread.
// The levels_15C analysis uses ebzSlope = 11 MeV/m, chosen empirically on the data.
double rmsAt(const char* f,double slope,double Ebeam=191.25){
  TFile fm(f); TTree* t=(TTree*)fm.Get("cbmsim");
  TClonesArray* mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
  Pars P; P.ebeam=Ebeam; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=true; P.ebzSlope=slope; P.ebzPivot=500;
  const double md=2.0141017781*931.49401;
  std::vector<double> v;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    for(int j=0;j<mc->GetEntriesFast();j++){
      auto* m=(FairMCTrack*)mc->At(j);
      if(m->GetPdgCode()!=1000010020||m->GetMotherId()!=-1) continue;
      double px=m->GetPx()*1000,py=m->GetPy()*1000,pz=m->GetPz()*1000;
      double p=std::sqrt(px*px+py*py+pz*pz), ke=std::sqrt(p*p+md*md)-md;
      double th=std::acos(pz/p)*180.0/M_PI, vz=m->GetStartZ()*10.0;
      double e=exOf(P,ebAt(P,vz,true),th,ke);
      if(std::isfinite(e)) v.push_back(e);
      break; } }
  if(v.size()<10) return -1;
  double s=0; for(double x:v) s+=x; double mu=s/v.size();
  double s2=0; for(double x:v) s2+=(x-mu)*(x-mu);
  return std::sqrt(s2/v.size());
}
void ebzscan(){
  const char* F[2]={"/home/yassid/Data/C15dd_sim/tc_0.0_sim.root",
                    "/home/yassid/Data/C15dd_sim/tc_6.68_sim.root"};
  double EX[2]={0.0,6.68};
  printf("  TRUTH-LEVEL Ex spread (rms, MeV) vs the Ebeam(vz) slope -- no detector, no reconstruction\n");
  printf("  slope:"); for(double s=0;s<=26.01;s+=2) printf("%7.0f",s); printf("   -> best\n");
  for(int i=0;i<2;i++){
    printf("  Ex %4.2f",EX[i]); double best=1e9,bs=0;
    for(double s=0;s<=26.01;s+=2){ double r=rmsAt(F[i],s); printf("%7.3f",r);
      if(r>0&&r<best){best=r;bs=s;} }
    printf("    %.0f MeV/m (rms %.3f)\n",bs,best); }
  printf("\n  the analysis uses ebzSlope = 11 MeV/m, set empirically on the DATA.\n");
}
