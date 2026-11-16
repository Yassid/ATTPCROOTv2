#include "binacc.h"
#include "binbuild.C"
// Double-scattering test.  An elastic deuteron starts ON the beam axis.  A deuteron produced by
// a SECOND scatter, tens of cm along the first deuteron's path, starts far off it.  Join the
// gated (d,d') events back to the upstream per-run fit trees, which carry vx,vy, and compare the
// transverse vertex radius of the background against the elastic peak.
void vxy(){
  setCM(10,5,60);
  // ---- upstream vx,vy, keyed by (run,event,trackID)
  TChain ch("kin"); ch.Add("/home/yassid/C15d_fit/*_kin_d.root");
  Int_t run,event,track; Double_t vx,vy,vzu,phiX;
  ch.SetBranchAddress("run",&run); ch.SetBranchAddress("event",&event);
  ch.SetBranchAddress("trackID",&track);
  ch.SetBranchAddress("vx",&vx); ch.SetBranchAddress("vy",&vy); ch.SetBranchAddress("vz",&vzu);
  ch.SetBranchAddress("phiXtr",&phiX);
  struct V{float x,y,phi;};
  std::unordered_map<long long,V> M; M.reserve(ch.GetEntries()*1.4);
  for(Long64_t i=0;i<ch.GetEntries();i++){ ch.GetEntry(i);
    M[((long long)run<<40)|((long long)event<<10)|(track&0x3FF)]={(float)vx,(float)vy,(float)phiX}; }
  printf("   upstream vertices loaded: %zu\n",M.size());
  // ---- beam axis = median of the elastic vertices
  Cfg& k=CFG[0]; Pars A=parsOf(k);
  TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2; float icf; int np,r2,e2,t2;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);t->SetBranchAddress("ic",&icf);t->SetBranchAddress("npulse",&np);
  t->SetBranchAddress("run",&r2);t->SetBranchAddress("event",&e2);t->SetBranchAddress("trackID",&t2);
  std::vector<float> gx,gy;
  struct S{const char* n; double lo,hi; std::vector<float> r;};
  S w[4]={{"g.s.   Ex -0.5..0.5",-0.5,0.5,{}},
          {"       Ex  1.0..3.0", 1.0,3.0,{}},
          {"bkg    Ex  3.0..5.0", 3.0,5.0,{}},
          {"bkg    Ex  7.0..9.0", 7.0,9.0,{}}};
  std::vector<std::pair<double,long long>> keep;  // (Ex, key)
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icf<A.icLo||icf>A.icHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
    double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>20) continue;
    long long key=((long long)r2<<40)|((long long)e2<<10)|(t2&0x3FF);
    auto it=M.find(key); if(it==M.end()) continue;
    if(e>-0.5&&e<0.5){ gx.push_back(it->second.x); gy.push_back(it->second.y); }
    keep.push_back({e,key}); }
  f->Close();
  auto med=[](std::vector<float>&v){ if(v.empty())return 0.f; std::nth_element(v.begin(),v.begin()+v.size()/2,v.end()); return v[v.size()/2]; };
  float cx=med(gx), cy=med(gy);
  printf("   beam axis from %zu elastic vertices: (x,y) = (%.1f, %.1f) mm\n\n",gx.size(),cx,cy);
  for(auto&p:keep){ auto&v=M[p.second];
    double r=std::hypot(v.x-cx,v.y-cy);
    for(auto&s:w) if(p.first>=s.lo&&p.first<s.hi) s.r.push_back(r); }
  printf("== ATTPCROOT theta_CM 10-20 deg : transverse vertex radius from the beam axis\n");
  printf("   %-22s %8s   %7s %7s %7s %7s   %8s %8s\n","window","N","med","68%","90%","95%","f(r>30)","f(r>60)");
  for(auto&s:w){ if(s.r.size()<20) continue; std::sort(s.r.begin(),s.r.end()); size_t N=s.r.size();
    auto q=[&](double fr){ size_t j=(size_t)(fr*N); if(j>=N) j=N-1; return s.r[j]; };
    size_t n30=0,n60=0; for(float v:s.r){ if(v>30) n30++; if(v>60) n60++; }
    printf("   %-22s %8zu   %7.1f %7.1f %7.1f %7.1f   %7.1f%% %7.1f%%\n",
           s.n,N,q(0.5),q(0.68),q(0.90),q(0.95),100.0*n30/N,100.0*n60/N); }
  printf("\n   Elastic deuterons start ON the beam axis; a second scatter happens tens of cm\n"
         "   along the first deuteron's path, so it would start far off it.\n");
}
