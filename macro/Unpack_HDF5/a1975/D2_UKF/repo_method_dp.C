/// Apply the notebook (draw_RCS.ipynb cell 7) method to OUR genfit data: counts in a vz slice, x (1/slice fraction),
/// divided by the repo efficiency table (efficiency_16Cdp_0m.txt), L = 316.98. Same bound selection as angdist_bound_dp.C.
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
void repo_method_dp(){
  double ang[6]={12.5,17.5,22.5,27.5,32.5,37.5}, effR[6]={0.828543,0.796656,0.860641,0.850787,0.824078,0.751247};
  double repoCnt[6]={72,70,21,25,32,27};
  double S[4][2]={{900,1000},{800,1000},{100,500},{500,800}}; const char*nm[4]={"900-1000","800-1000","100-500","500-800"};
  double C[4][6]={{0}};
  TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) continue;
    double ex,tc; if(!kine2b(184.25-0.011030*(vz-500),th*TMath::DegToRad(),ke,ex,tc)) continue; if(ex<-1.5||ex>0.729) continue;
    int b=(int)((tc-10)/5); if(b<0||b>5) continue; for(int s=0;s<4;++s) if(vz>=S[s][0]&&vz<S[s][1]) C[s][b]++; }
  printf("counts per theta bin by vz slice (our genfit data) vs the notebook's Spyral counts\n  thcm  ");
  for(int s=0;s<4;++s) printf("  %9s",nm[s]); printf("   | notebook 0.9-1.0 m\n");
  for(int b=0;b<6;++b){ printf("  %4.1f  ",ang[b]); for(int s=0;s<4;++s) printf("  %9.0f",C[s][b]); printf("   | %5.0f\n",repoCnt[b]); }
  printf("\nNOTEBOOK METHOD: sigma = counts/(L*dOmega) x (1/slice fraction) / eff_repo, L = 316.98\n  thcm    ours 900-1000 x10   notebook (cell 7)   ours full 0-915 (sim acceptance)\n");
  double full[6]={1.637,1.285,1.686,2.994,6.666,3.441};
  for(int b=0;b<6;++b){ double dO=2*TMath::Pi()*(cos((ang[b]-2.5)*TMath::DegToRad())-cos((ang[b]+2.5)*TMath::DegToRad()));
    double ours=C[0][b]/(316.98*dO)*10/effR[b], rep=repoCnt[b]/(316.98*dO)*10/effR[b];
    printf("  %4.1f      %6.2f +- %4.2f        %6.2f              %6.2f\n",ang[b],ours,ours/sqrt(std::max(1.,C[0][b])),rep,full[b]); }
}
