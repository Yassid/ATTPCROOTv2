/// @file radius_bias_check.C
/// @brief Localise the momentum bias: compare, at fixed |p|, three radii:
///   R0    = analytic vertex radius p/(0.3 B)
///   Rtru  = circle fit to the TRUTH MC points (x,y)  -> energy-loss spiral avg
///   Rreco = PRA GetGeoRadius (reco hits)             -> + PSA/clustering effects
/// Rtru/R0-1 is the energy-loss contribution; Rreco/Rtru-1 is the reco contribution.
/// Run: root -b -q 'radius_bias_check.C("data/attpcsim.root","data/scan_p0.900.root",0.900)'
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

// algebraic (Kasa) circle fit -> radius [mm]
double circR(const std::vector<double>&x,const std::vector<double>&y){
   int n=x.size(); if(n<4)return 0;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxz=0,Syz=0,Sz=0;
   for(int i=0;i<n;i++){ double z=x[i]*x[i]+y[i]*y[i];
      Sx+=x[i];Sy+=y[i];Sxx+=x[i]*x[i];Syy+=y[i]*y[i];Sxy+=x[i]*y[i];Sxz+=x[i]*z;Syz+=y[i]*z;Sz+=z; }
   // solve [Sxx Sxy Sx; Sxy Syy Sy; Sx Sy n] [D E F]' = -[Sxz Syz Sz]'
   double a[3][3]={{Sxx,Sxy,Sx},{Sxy,Syy,Sy},{Sx,Sy,(double)n}};
   double b[3]={-Sxz,-Syz,-Sz};
   // Gaussian elimination
   for(int c=0;c<3;c++){ int p=c; for(int r=c+1;r<3;r++) if(fabs(a[r][c])>fabs(a[p][c]))p=r;
      std::swap(a[c],a[p]); std::swap(b[c],b[p]); if(fabs(a[c][c])<1e-12)return 0;
      for(int r=0;r<3;r++) if(r!=c){ double f=a[r][c]/a[c][c]; for(int k=0;k<3;k++)a[r][k]-=f*a[c][k]; b[r]-=f*b[c]; } }
   double D=b[0]/a[0][0],E=b[1]/a[1][1],F=b[2]/a[2][2];
   double disc=D*D/4+E*E/4-F; return disc>0?std::sqrt(disc):0;
}

void radius_bias_check(TString simFile="data/attpcsim.root", TString recoFile="data/scan_p0.900.root", double pGeV=0.900)
{
   gSystem->Load("libAtReconstruction.so");
   const double B=4.0, R0=pGeV/(0.2997925*B)*1000.0; // mm
   // ---- truth radii from MC points ----
   TFile fS(simFile); auto*tS=(TTree*)fS.Get("cbmsim");
   auto*mcPts=new TClonesArray("AtMCPoint"); auto*mcTrks=new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint",&mcPts); tS->SetBranchAddress("MCTrack",&mcTrks);
   std::vector<double> Rtru;
   for(Long64_t e=0;e<tS->GetEntries();e++){ tS->GetEntry(e);
      std::map<int,std::vector<double>> xs,ys;
      for(int k=0;k<mcPts->GetEntries();k++){ auto*mp=(AtMCPoint*)mcPts->At(k);
         int tid=mp->GetTrackID(); auto*mt=(tid>=0&&tid<mcTrks->GetEntries())?(AtMCTrack*)mcTrks->At(tid):nullptr;
         if(!mt||abs(mt->GetPdgCode())!=211)continue;                 // primary pions only
         xs[tid].push_back(mp->GetX()*10); ys[tid].push_back(mp->GetY()*10); }
      for(auto&kv:xs){ double r=circR(kv.second,ys[kv.first]); if(r>50&&r<5000)Rtru.push_back(r); }
   }
   // ---- reco radii ----
   TFile fD(recoFile); auto*tD=(TTree*)fD.Get("cbmsim");
   auto*pat=new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent",&pat);
   std::vector<double> Rreco, RrecoHit;
   for(Long64_t e=0;e<tD->GetEntries();e++){ tD->GetEntry(e); if(pat->GetEntries()==0)continue;
      auto*pe=(AtPatternEvent*)pat->At(0);
      for(auto&tr:pe->GetTrackCand()){ double R=tr.GetGeoRadius(); if(R>50&&R<5000)Rreco.push_back(R);
         std::vector<double> x,y; for(auto&h:tr.GetHitArray()){x.push_back(h->GetPosition().X());y.push_back(h->GetPosition().Y());}
         double rh=circR(x,y); if(rh>50&&rh<5000)RrecoHit.push_back(rh); }
   }
   double rt=med(Rtru), rr=med(Rreco), rrh=med(RrecoHit);
   printf("\n==== radius bias @ |p|=%.3f GeV  (R0=%.1f mm) ====\n",pGeV,R0);
   printf("  R_truth(MC circle)   = %.1f mm   (%+.1f%% vs R0)   [energy-loss spiral]\n",rt,100*(rt-R0)/R0);
   printf("  R_reco (GetGeoRadius)= %.1f mm   (%+.1f%% vs R0,  %+.1f%% vs Rtru)\n",rr,100*(rr-R0)/R0,rt>0?100*(rr-rt)/rt:0);
   printf("  R_reco (refit hits)  = %.1f mm   (%+.1f%% vs R0,  %+.1f%% vs Rtru)\n",rrh,100*(rrh-R0)/R0,rt>0?100*(rrh-rt)/rt:0);
   printf("  => energy-loss part %+.1f%% ; reco/PSA part %+.1f%%\n",100*(rt-R0)/R0, rt>0?100*(rr-rt)/rt:0);
}
