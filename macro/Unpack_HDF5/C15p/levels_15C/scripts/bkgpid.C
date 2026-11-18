#include "binacc.h"
#include "binbuild.C"
// IS THE BACKGROUND PROTON OR DEUTERON?
// The PID plane is ln(sqrt(dEdx) * Brho).  At the SAME Brho a proton and a deuteron carry the same
// momentum, so the deuteron is slower and loses MORE energy: deuterons sit ABOVE protons.
// Measure the locus of the background against the elastic (certain deuterons) and against the
// proton band, on the clean window and the loosened cut.
void bkgpid(){
  TFile* fp=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root");
  TTree* pt=(TTree*)fp->Get("pts");
  float sq,br; int rn,ev,tid;
  pt->SetBranchAddress("sqrtdedx",&sq); pt->SetBranchAddress("brho",&br);
  pt->SetBranchAddress("run",&rn); pt->SetBranchAddress("event",&ev); pt->SetBranchAddress("trackID",&tid);
  std::unordered_map<long long,float> M; M.reserve(pt->GetEntries()*1.4);
  // the proton / deuteron loci straight from the plane, before any gate
  TH1D* hall=new TH1D("hall","",400,1.0,4.0); hall->SetDirectory(nullptr);
  for(Long64_t i=0;i<pt->GetEntries();i++){ pt->GetEntry(i);
    if(sq<=0||br<=0) continue;
    double L=std::log(sq*br);
    hall->Fill(L);
    M[((long long)rn<<40)|((long long)ev<<10)|(tid&0x3FF)]=(float)L; }
  fp->Close();
  printf("\n================ the PID plane itself (ungated) ================\n");
  printf("   ln(sqrt(dEdx)*Brho): the two loci, found as the two largest maxima\n");
  std::vector<std::pair<double,double>> pk;
  for(int i=3;i<hall->GetNbinsX()-2;i++){
    double y=hall->GetBinContent(i);
    if(y>hall->GetBinContent(i-1)&&y>hall->GetBinContent(i+1)&&
       y>hall->GetBinContent(i-2)&&y>hall->GetBinContent(i+2)&&y>hall->GetEntries()*0.002)
      pk.push_back({y,hall->GetBinCenter(i)}); }
  std::sort(pk.rbegin(),pk.rend());
  for(size_t i=0;i<std::min<size_t>(4,pk.size());i++)
    printf("     locus %zu at ln = %.4f   (height %.0f)\n",i+1,pk[i].second,pk[i].first);
  double sep = (pk.size()>1)? std::fabs(pk[0].second-pk[1].second) : 0;
  printf("   separation between the two strongest: %.4f\n",sep);

  for(int ic=0;ic<1;ic++){
    Cfg k=CFG[ic]; k.vzHi=900; k.chi2=5.0;          // the loosened cut
    Pars A=parsOf(k); A.vzLo=300;
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0; int r2,e2,t2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("run",&r2);t->SetBranchAddress("event",&e2);t->SetBranchAddress("trackID",&t2);
    struct W{const char* n; double lo,hi; long N; double s,s2;};
    W ws[7]={{"elastic   |Ex|<0.5",-0.5,0.5,0,0,0},
             {"Ex 2-3            ", 2.0,3.0,0,0,0},
             {"Ex 3-4            ", 3.0,4.0,0,0,0},
             {"Ex 4-5            ", 4.0,5.0,0,0,0},
             {"Ex 5-7            ", 5.0,7.0,0,0,0},
             {"Ex 7-9            ", 7.0,9.0,0,0,0},
             {"Ex 9-12           ", 9.0,12.0,0,0,0}};
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
      auto it=M.find(((long long)r2<<40)|((long long)e2<<10)|(t2&0x3FF));
      if(it==M.end()) continue;
      for(auto&wv:ws) if(e>=wv.lo&&e<wv.hi){ wv.N++; wv.s+=it->second; wv.s2+=it->second*it->second; } }
    f->Close();
    printf("\n================ %s : where each Ex window sits in the PID plane ================\n",k.nm);
    printf("   %-20s %9s %12s %9s %14s\n","window","N","<ln>","rms","shift vs elastic");
    double ref=0;
    for(auto&wv:ws){ if(!wv.N) continue;
      double m=wv.s/wv.N, r=std::sqrt(std::max(0.0,wv.s2/wv.N-m*m));
      if(ref==0) ref=m;
      printf("   %-20s %9ld %12.4f %9.4f %+14.4f\n",wv.n,wv.N,m,r,m-ref); }
    printf("\n   A proton contamination would pull the mean DOWN (protons sit below deuterons).\n");
    printf("   Scale: the two loci are %.4f apart.\n",sep);
  }
}
