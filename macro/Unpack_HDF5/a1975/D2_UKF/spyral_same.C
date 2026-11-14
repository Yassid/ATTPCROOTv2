// Spyral-1.0 InterpSolver (d,p) data through EXACTLY the genfit bound-state method of angdist_bound_dp.C
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz){ return 184.25-0.011030*(vz-500); }
struct Ev{double ex,tc,vz;};
const int NT=7; double TB[NT+1]={10,15,20,25,30,35,40,45};
int tb(double t){ for(int i=0;i<NT;++i) if(t>=TB[i]&&t<TB[i+1]) return i; return -1; }
void peak(const char*nm,std::vector<Ev>&v){ // robust core of the bound peak, vz 200-400, thcm 15-40
  std::vector<double> x; for(auto&e:v) if(e.vz>=200&&e.vz<400&&e.tc>=15&&e.tc<40&&e.ex>-2.5&&e.ex<1.5) x.push_back(e.ex);
  double m=0,s=0.4; for(int it=0;it<10;++it){ std::vector<double> y; for(double a:x) if(fabs(a-m)<2.5*s) y.push_back(a); std::sort(y.begin(),y.end());
    m=y[y.size()/2]; s=0.5*(y[(int)(0.84*y.size())]-y[(int)(0.16*y.size())]); }
  printf("%-8s bound-peak core (vz 200-400, thcm 15-40): centroid %+.3f MeV  sigma %.3f  n %zu\n",nm,m,s,x.size()); }
void spyral_same(double shiftSpy=0){
  std::vector<Ev> G,S;
  { TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree*t=(TTree*)f.Get("pk"); float ke,th,vz,c2;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vertexz",&vz);t->SetBranchAddress("chi2ndf",&c2);
    for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); if(!(c2<=13.5)||!(th>=89&&th<=161)||!(ke>=1.2)) continue; double ex,tc; if(!kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc)) continue; G.push_back({ex,tc,vz});} }
  long nS=0;
  { void*d=gSystem->OpenDirectory("/home/yassid/C16_dp/C16_dp/InterpSolver_root"); const char*fn;
    while((fn=gSystem->GetDirEntry(d))){ TString s(fn); if(!s.EndsWith("_1H.root")) continue;
      TFile f("/home/yassid/C16_dp/C16_dp/InterpSolver_root/"+s); TTree*t=(TTree*)f.Get("parquettree"); double pol,br,z,rc;
      t->SetBranchAddress("polar",&pol);t->SetBranchAddress("brho",&br);t->SetBranchAddress("vertex_z",&z);t->SetBranchAddress("redchisq",&rc);
      for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i); nS++;
        double p=br*2.99792458e2, ke=sqrt(p*p+M3*M3)-M3, th=pol*TMath::RadToDeg(), vz=z*1000;
        if(!(th>=89&&th<=161)||!(ke>=1.2)) continue; double ex,tc; if(!kine2b(Eb(vz),th*TMath::DegToRad(),ke,ex,tc)) continue; S.push_back({ex+shiftSpy,tc,vz}); } }
    gSystem->FreeDirectory(d); }
  printf("Spyral tracks read %ld; after theta_lab 89-161, KE>1.2: Spyral %zu  genfit %zu\n",nS,S.size(),G.size());
  peak("genfit",G); peak("Spyral",S);
  // counts per theta bin by vz slice (bound window + sideband subtracted)
  double VZ[]={0,100,200,300,400,500,600,700,800,900,1000}; const int NV=10;
  auto cnt=[&](std::vector<Ev>&v,int b,double z0,double z1){ double y=0,bg=0; for(auto&e:v){ if(e.vz<z0||e.vz>=z1||tb(e.tc)!=b) continue;
      if(e.ex>=-1.5&&e.ex<=0.729) y++; else if(e.ex>=-4&&e.ex<-1.5) bg++; } return y-bg*2.229/2.5; };
  printf("\nnet bound counts, genfit / Spyral, per 100 mm vz slice\n thcm ");
  for(int v=0;v<NV;++v) printf(" %4.0f-%-4.0f ",VZ[v],VZ[v+1]); printf("\n");
  for(int b=0;b<NT;++b){ printf("%2.0f-%2.0f",TB[b],TB[b+1]); for(int v=0;v<NV;++v) printf("  %4.0f/%-4.0f",cnt(G,b,VZ[v],VZ[v+1]),cnt(S,b,VZ[v],VZ[v+1])); printf("\n"); }
  // dsdo in vz 200-400 with the SAME sim eps (dp_bound_angdist_vz200-400.csv) and L_W
  double eps[NT]={0.8842,0.8954,0.7763,0.9023,0.7021,0.7695,0.7116}, LW=316.98*200/971.73;
  printf("\ndsigma/dOmega [mb/sr], vz 200-400, L_W %.2f, same g.s.-sim eps for both\n thcm   Y_gf  Y_spy   genfit          Spyral          Spy/gf\n",LW);
  FILE*fo=fopen("spyral_same_vz200-400.csv","w"); fprintf(fo,"thcm_lo,thcm_hi,Ygf,Yspy,dsdo_gf,err_gf,dsdo_spy,err_spy\n");
  for(int b=0;b<NT;++b){ double dO=2*TMath::Pi()*(cos(TB[b]*TMath::DegToRad())-cos(TB[b+1]*TMath::DegToRad()));
    double yg=cnt(G,b,200,400), ys=cnt(S,b,200,400), k=1/(eps[b]*dO*LW);
    printf("%2.0f-%2.0f  %5.0f  %5.0f   %6.2f(%4.2f)    %6.2f(%4.2f)    %5.2f\n",TB[b],TB[b+1],yg,ys,yg*k,sqrt(std::max(1.,yg))*k,ys*k,sqrt(std::max(1.,ys))*k,yg>0?ys/yg:0);
    fprintf(fo,"%.0f,%.0f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4f\n",TB[b],TB[b+1],yg,ys,yg*k,sqrt(std::max(1.,yg))*k,ys*k,sqrt(std::max(1.,ys))*k); }
  fclose(fo);
}
