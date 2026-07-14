/// @file ring_centroid_radius.C
/// @brief Is the DLC momentum bias in the digi (charge centroid shifted) or the
///        reco (pad-quantized hits)? Per PRA track, fit a circle to (a) the raw
///        pad-centre hits and (b) the charge-weighted per-ring centroids, and
///        compare the radius. If (b) >> (a) and near truth, charge-centroiding
///        fixes it and should be pushed into pattern recognition, not just the fit.
/// Run: root -b -q 'ring_centroid_radius.C("data/dlc_on.root","DLC on")'
double medd(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
// algebraic (Kasa) circle radius [mm]
double kasaR(const std::vector<double>&x,const std::vector<double>&y){
   int n=x.size(); if(n<4)return 0;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxz=0,Syz=0,Sz=0;
   for(int i=0;i<n;i++){ double z=x[i]*x[i]+y[i]*y[i];
      Sx+=x[i];Sy+=y[i];Sxx+=x[i]*x[i];Syy+=y[i]*y[i];Sxy+=x[i]*y[i];Sxz+=x[i]*z;Syz+=y[i]*z;Sz+=z; }
   double a[3][3]={{Sxx,Sxy,Sx},{Sxy,Syy,Sy},{Sx,Sy,(double)n}}, b[3]={-Sxz,-Syz,-Sz};
   for(int c=0;c<3;c++){ int p=c; for(int r=c+1;r<3;r++) if(fabs(a[r][c])>fabs(a[p][c]))p=r;
      std::swap(a[c],a[p]); std::swap(b[c],b[p]); if(fabs(a[c][c])<1e-12)return 0;
      for(int r=0;r<3;r++) if(r!=c){ double f=a[r][c]/a[c][c]; for(int k=0;k<3;k++)a[r][k]-=f*a[c][k]; b[r]-=f*b[c]; } }
   double D=b[0]/a[0][0],E=b[1]/a[1][1],F=b[2]/a[2][2]; double d=D*D/4+E*E/4-F; return d>0?sqrt(d):0;
}

void ring_centroid_radius(TString recoFile="data/dlc_on.root", TString tag="DLC on", int nP=256)
{
   gSystem->Load("libAtReconstruction.so");
   TFile f(recoFile); auto*t=(TTree*)f.Get("cbmsim");
   TClonesArray*pat=nullptr; t->SetBranchAddress("AtPatternEvent",&pat);
   std::vector<double> Rraw, Rcen;
   for(Long64_t e=0;e<t->GetEntries();e++){ t->GetEntry(e); if(!pat->GetEntries())continue;
      for(auto&tr:((AtPatternEvent*)pat->At(0))->GetTrackCand()){
         std::vector<double> xr,yr;
         std::map<int,std::vector<int>> byRing; auto&h=tr.GetHitArray();
         for(int i=0;i<(int)h.size();i++){ const auto&p=h[i]->GetPosition(); xr.push_back(p.X()); yr.push_back(p.Y());
            byRing[h[i]->GetPadNum()/nP].push_back(i); }
         double rr=kasaR(xr,yr); if(rr>50&&rr<5000)Rraw.push_back(rr);
         // charge-weighted centroid per ring
         std::vector<double> xc,yc;
         for(auto&kv:byRing){ double sx=0,sy=0,sq=0; for(int i:kv.second){ double q=h[i]->GetCharge();
            sx+=h[i]->GetPosition().X()*q; sy+=h[i]->GetPosition().Y()*q; sq+=q; }
            if(sq>0){ xc.push_back(sx/sq); yc.push_back(sy/sq); } }
         double rc=kasaR(xc,yc); if(rc>50&&rc<5000)Rcen.push_back(rc);
      }
   }
   printf("  %-10s  R_raw(pad-centre)=%.1f mm   R_cen(charge-weighted rings)=%.1f mm   [truth R0(375)=312, in-gas~299]\n",
          tag.Data(), medd(Rraw), medd(Rcen));
}
