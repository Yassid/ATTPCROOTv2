#include "/home/yassid/Data/C15d_statefit/excm.h"
// Why does the acceptance denominator keep only ~half the generated events?
void denom(){
  const double md=2.0141017781*931.49401;
  Pars P; P.ebeam=191.25; P.mBeam=15.0105993; P.mTarg=2.014102; P.mEject=2.014102; P.mResid=15.0105993;
  P.kcOn=false; P.ebzOn=true; P.ebzSlope=11.0; P.ebzPivot=500;
  TFile fm("/home/yassid/Data/C15dd_sim/ex00.00_sim.root");
  TTree* t=(TTree*)fm.Get("cbmsim");
  TClonesArray* mc=nullptr; t->SetBranchAddress("MCTrack",&mc);
  long tot=0,noD=0,badVz=0,badCm=0,ok=0;
  TH1D hv("hv","",60,-200,1400), hc("hc","",90,0,180);
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i); tot++;
    bool found=false; double vz=0,tcm=-1;
    for(int j=0;j<mc->GetEntriesFast();j++){
      auto* m=(FairMCTrack*)mc->At(j);
      if(m->GetPdgCode()!=1000010020||m->GetMotherId()!=-1) continue;
      double px=m->GetPx()*1000,py=m->GetPy()*1000,pz=m->GetPz()*1000;
      double p=std::sqrt(px*px+py*py+pz*pz), ke=std::sqrt(p*p+md*md)-md;
      vz=m->GetStartZ()*10.0;
      tcm=thcmProj(P,ebAt(P,vz,true),std::acos(pz/p)*180.0/M_PI,ke);
      found=true; break; }
    if(!found){noD++;continue;}
    hv.Fill(vz); if(std::isfinite(tcm)) hc.Fill(tcm);
    if(!std::isfinite(tcm)||tcm<0){badCm++;continue;}
    if(vz<0||vz>1100){badVz++;continue;}
    ok++; }
  printf("  total %ld | no primary deuteron %ld | theta_CM invalid %ld | vz outside [0,1100] %ld | kept %ld (%.1f%%)\n",
         tot,noD,badCm,badVz,ok,100.0*ok/tot);
  printf("  vz   : "); for(double x=-100;x<=1300;x+=200) printf(" %5.0f:%5.0f",x,hv.GetBinContent(hv.FindBin(x))); printf("\n");
  printf("  thcm : "); for(double x=5;x<=175;x+=20) printf(" %5.0f:%5.0f",x,hc.GetBinContent(hc.FindBin(x))); printf("\n");
}
