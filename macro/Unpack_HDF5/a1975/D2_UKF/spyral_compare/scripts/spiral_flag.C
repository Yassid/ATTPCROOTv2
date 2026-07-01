// Per-event adaptive trigger inputs: largest-track size in the default and loose
// PRA recos, plus the turn-count of the loose largest track (multi-turn spiral
// detector). The adaptive merge (Python) flags an event to use the LOOSE clustering
// when loose recovered substantially more AND the track is a genuine multi-turn
// spiral (so messy proton+beam over-merges are NOT flagged).
//
//   root -l -b -q 'spiral_flag.C("run_0016")'  -> spyral_compare/spiral_flag.csv
#include <vector>
#include <algorithm>

static double track_turns(AtTrack &tr) {
   auto &hits = tr.GetHitArray();
   int n = hits.size();
   if (n < 10) return 0;
   std::vector<double> x(n), y(n), z(n);
   for (int i = 0; i < n; ++i) { auto p = hits[i]->GetPosition(); x[i]=p.X(); y[i]=p.Y(); z[i]=p.Z(); }
   // Kasa circle center
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxz=0,Syz=0,Szz=0;
   for (int i=0;i<n;++i){ double zz=x[i]*x[i]+y[i]*y[i];
      Sx+=x[i];Sy+=y[i];Sxx+=x[i]*x[i];Syy+=y[i]*y[i];Sxy+=x[i]*y[i];Sxz+=x[i]*zz;Syz+=y[i]*zz; }
   double A[2][2]={{Sxx,Sxy},{Sxy,Syy}}, b[2]={0.5*Sxz,0.5*Syz};
   // subtract means for stability
   double mx=Sx/n,my=Sy/n; A[0][0]=Sxx-n*mx*mx; A[0][1]=A[1][0]=Sxy-n*mx*my; A[1][1]=Syy-n*my*my;
   double bz0=0,bz1=0,mz=0; for(int i=0;i<n;++i){double zz=x[i]*x[i]+y[i]*y[i];mz+=zz;} mz/=n;
   for(int i=0;i<n;++i){double zz=x[i]*x[i]+y[i]*y[i]-mz;bz0+=(x[i]-mx)*zz;bz1+=(y[i]-my)*zz;}
   double det=A[0][0]*A[1][1]-A[0][1]*A[1][0]; if(std::fabs(det)<1e-9)return 0;
   double cx=mx+0.5*( A[1][1]*bz0-A[0][1]*bz1)/det;
   double cy=my+0.5*(-A[1][0]*bz0+A[0][0]*bz1)/det;
   // order by z, unwrap angle
   std::vector<int> idx(n); for(int i=0;i<n;++i)idx[i]=i;
   std::sort(idx.begin(),idx.end(),[&](int a,int b){return z[a]<z[b];});
   double prev=0,cum=0,amin=1e9,amax=-1e9; bool first=true;
   for(int k=0;k<n;++k){int i=idx[k];double a=std::atan2(y[i]-cy,x[i]-cx);
      if(first){prev=a;cum=a;first=false;}
      else{double d=a-prev; while(d>M_PI)d-=2*M_PI; while(d<-M_PI)d+=2*M_PI; cum+=d; prev=a;}
      amin=std::min(amin,cum); amax=std::max(amax,cum);}
   return (amax-amin)/(2*M_PI);
}

static int maxTrack(TClonesArray* pat, double* turnsOut=nullptr){
   if(!pat||pat->GetEntriesFast()==0){if(turnsOut)*turnsOut=0;return 0;}
   auto pe=(AtPatternEvent*)pat->At(0); int mx=0; AtTrack* big=nullptr;
   for(auto& tr:pe->GetTrackCand()){int s=tr.GetHitArray().size(); if(s>mx){mx=s;big=&tr;}}
   if(turnsOut)*turnsOut= big? track_turns(*big):0;
   return mx;
}

void spiral_flag(TString run="run_0016", TString ioDir="/mnt/f/a1975/reco_d2/"){
   TFile* fd=TFile::Open(ioDir+run+"_reco.root");
   TFile* fl=TFile::Open(ioDir+run+"loose_reco.root");
   TTree* td=(TTree*)fd->Get("cbmsim"); TTree* tl=(TTree*)fl->Get("cbmsim");
   TClonesArray *pd=0,*pl=0; td->SetBranchAddress("AtPatternEvent",&pd); tl->SetBranchAddress("AtPatternEvent",&pl);
   Long64_t N=std::min(td->GetEntries(),tl->GetEntries());
   FILE* fp=fopen("spyral_compare/spiral_flag.csv","w");
   fprintf(fp,"event,def_max,loose_max,loose_turns\n");
   for(Long64_t i=0;i<N;++i){ td->GetEntry(i); tl->GetEntry(i);
      double turns=0; int lm=maxTrack(pl,&turns); int dm=maxTrack(pd);
      fprintf(fp,"%lld,%d,%d,%.2f\n",i,dm,lm,turns);
      if(i%5000==0)printf("  %lld/%lld\n",i,N);
   }
   fclose(fp); printf("DONE -> spyral_compare/spiral_flag.csv (%lld events)\n",N);
}
