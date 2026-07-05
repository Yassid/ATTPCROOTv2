/// @file pra_analyze.C
/// @brief Evaluate PRA-setting variants for the 200 MeV/c pi over-segmentation problem.
///        For each config file: PRA tracks/event, momentum-rigidity shape (skew, exKurt,
///        IQR, median bias), and the fraction of tracks whose radius is within 30% of the
///        truth (167 mm). Goal: fewer tracks (toward 2), lower skew/kurt (toward Gaussian).
/// Run: root -b -q pra_analyze.C
double skw(std::vector<double>v){if(v.size()<3)return 0;double m=0;for(double x:v)m+=x;m/=v.size();double s2=0,s3=0;for(double x:v){double d=x-m;s2+=d*d;s3+=d*d*d;}s2/=v.size();s3/=v.size();return s2>0?s3/pow(s2,1.5):0;}
double krt(std::vector<double>v){if(v.size()<4)return 0;double m=0;for(double x:v)m+=x;m/=v.size();double s2=0,s4=0;for(double x:v){double d=x-m;s2+=d*d;s4+=d*d*d*d;}s2/=v.size();s4/=v.size();return s2>0?s4/(s2*s2)-3:0;}
double iqrp(std::vector<double>v){if(v.size()<4)return 0;std::sort(v.begin(),v.end());return(v[3*v.size()/4]-v[v.size()/4])/1.349;}
double medp(std::vector<double>v){if(v.empty())return 0;std::sort(v.begin(),v.end());return v[v.size()/2];}

void eval(const char*tag,TString file,double B=4.0,double p0=200,double Rtruth=166.8){
   TFile f(file); TTree*t=(TTree*)f.Get("cbmsim");
   if(!t){printf("  %-8s MISSING\n",tag);return;}
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&pat);
   long nev=0,ntrk=0; std::vector<double> dp; int good=0;
   for(Long64_t e=0;e<t->GetEntries();++e){t->GetEntry(e); if(!pat->GetEntries())continue; nev++;
      auto&tc=((AtPatternEvent*)pat->At(0))->GetTrackCand(); ntrk+=tc.size();
      for(auto&tr:tc){double R=tr.GetGeoRadius(); if(R>0&&R<1e5){dp.push_back(100*(0.299792458*B*R-p0)/p0);
         if(fabs(R-Rtruth)/Rtruth<0.3)good++;}} }
   printf("  %-8s tracks/evt=%.2f  n=%zu  skew=%+6.1f kurt=%+8.0f IQR=%5.1f%% med=%+5.1f%%  goodR=%.0f%%\n",
      tag,nev?(double)ntrk/nev:0,dp.size(),skw(dp),krt(dp),iqrp(dp),medp(dp),dp.empty()?0:100.0*good/dp.size());
}

void pra_analyze(){
   gSystem->Load("libAtReconstruction.so");
   TString D="/mnt/f/puma_sweep/";
   printf("\n===== PRA-setting scan @ 200 MeV/c (truth: 2 tracks, R=167mm; want low skew/kurt) =====\n");
   eval("baseline",D+"output_digi_pi200_primary.root"); // merge off, no minHits
   eval("merge",   D+"output_digi_pra_merge.root");
   eval("min15",   D+"output_digi_pra_min15.root");
   eval("min15+mrg",D+"output_digi_pra_min15m.root");
   eval("min25",   D+"output_digi_pra_min25.root");
   printf("=================================================================================\n\n");
}
