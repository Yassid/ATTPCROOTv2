/// @file gf_hyp_sep_dp.C
/// @brief d/t misidentification per circle-Brho slice for the genfit two-hypothesis test (2026-11-01), on the d + t
///        sims: (1) dE/dx + circle Brho, (2) dE/dx + genfit d-hypothesis Brho, (3) chi2 difference D alone,
///        (4) (1) AND (3). Inputs: pid_sbpid_dp.C sim trees + gf_hyp_read_dp.C outputs in G (rows aligned).
///        Run with G = /mnt/f/a1975/gf_hyp/ (measSigma 4) or /mnt/f/a1975/gf_hyp/s06/ (0.6 mm). Fits: gf_hyp_run_dp.sh.
///   root -l -b -q 'gf_hyp_sep_dp.C("/mnt/f/a1975/gf_hyp/s06/")'
// per circle-Brho (bA) slice, d sim vs t sim: misID for
//  (1) dE/dx + circle Brho: cut on ln(hs*bA) at the width-weighted midpoint      [current adopted]
//  (2) dE/dx + genfit d-hypothesis Brho ln(hs*bx_d)
//  (3) genfit chi2 difference alone, D = (chi2_t - chi2_d)/ndf_d, cut D < Dc (best balanced Dc scanned)
//  (4) (1) AND (3): triton if ln-mass above the cut AND D < Dc
void gf_hyp_sep_dp(TString G="/mnt/f/a1975/gf_hyp/"){
  TChain d("w"),t("w"),dg("g"),tg("g");
  for(int s=5001;s<=5004;s++){ d.Add(Form("/mnt/f/a1975_C16_dd_sim/pv_el_s%d.root",s)); dg.Add(Form("%shyp_el_s%d.root",G.Data(),s)); }
  for(const char*s:{"gs","ex2"}){ t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root",s)); tg.Add(Form("%shyp_%s.root",G.Data(),s)); }
  d.AddFriend(&dg,"g"); t.AddFriend(&tg,"g");
  TTree*T[2]={&d,&t};
  TString sel="va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&bA>0&&g.ok_d&&g.ok_t";
  double bl[]={0.20,0.25,0.30,0.35,0.45,0.60,0.75};
  printf("bA slice     N_d  N_t | (1) circle lnm: t>d d>t | (2) genfit-d lnm: t>d d>t | (3) D only: Dc  t>d d>t | (4) (1)&(3): t>d d>t\n");
  for(int k=0;k<6;k++){
    std::vector<std::array<double,3>> v[2];
    for(int s=0;s<2;s++){ long n=T[s]->Draw("log(hs*bA):log(hs*g.bx_d):(g.chi2_t-g.chi2_d)/g.ndf_d",sel+Form("&&bA>%g&&bA<%g",bl[k],bl[k+1]),"goff");
      for(long i=0;i<n;i++) v[s].push_back({T[s]->GetV1()[i],T[s]->GetV2()[i],T[s]->GetV3()[i]}); }
    if(v[0].size()<15||v[1].size()<15){printf("%.2f-%.2f %4zu %4zu\n",bl[k],bl[k+1],v[0].size(),v[1].size());continue;}
    auto cutOf=[&](int c){ std::vector<double>a,b; for(auto&x:v[0])a.push_back(x[c]); for(auto&x:v[1])b.push_back(x[c]); std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
      auto q=[](std::vector<double>&z,double f){return z[std::min(z.size()-1,(size_t)(f*z.size()))];};
      double md=q(a,.5),mt=q(b,.5),hd=.5*(q(a,.84)-q(a,.16)),ht=.5*(q(b,.84)-q(b,.16)); return md+(mt-md)*hd/(hd+ht);};
    double c0=cutOf(0), c1=cutOf(1);
    auto frac=[&](int s,std::function<bool(const std::array<double,3>&)> isT){ long n=0; for(auto&x:v[s]) n+=isT(x); return 100.*n/v[s].size(); };
    // (3) scan Dc for the minimum of max(t>d, d>t)
    double bestDc=0,best=1e9; for(double Dc:{-1.0,-0.3,-0.1,-0.05,-0.02,-0.01,-0.005,-0.003,-0.002,-0.001,-0.0005,-0.0002,0.0}){
      double a=100-frac(1,[&](auto&x){return x[2]<Dc;}), b=frac(0,[&](auto&x){return x[2]<Dc;}); if(std::max(a,b)<best){best=std::max(a,b);bestDc=Dc;} }
    printf("%.2f-%.2f %4zu %4zu |  %5.1f%% %5.2f%%        |  %5.1f%% %5.2f%%           | %+.4f %5.1f%% %5.1f%%   | %5.1f%% %5.2f%%\n",bl[k],bl[k+1],v[0].size(),v[1].size(),
      100-frac(1,[&](auto&x){return x[0]>c0;}), frac(0,[&](auto&x){return x[0]>c0;}),
      100-frac(1,[&](auto&x){return x[1]>c1;}), frac(0,[&](auto&x){return x[1]>c1;}),
      bestDc, 100-frac(1,[&](auto&x){return x[2]<bestDc;}), frac(0,[&](auto&x){return x[2]<bestDc;}),
      100-frac(1,[&](auto&x){return x[0]>c0&&x[2]<bestDc;}), frac(0,[&](auto&x){return x[0]>c0&&x[2]<bestDc;}));
  }
}
