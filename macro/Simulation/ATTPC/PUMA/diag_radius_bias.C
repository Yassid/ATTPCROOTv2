/// @file diag_radius_bias.C
/// @brief Localize the radius UNDER-estimate on tight annular arcs. For truth-matched
///        pion tracks, compare three radius estimates against the truth:
///          R_PRA    — AtTrack::GetGeoRadius (cluster-based, what the pipeline uses)
///          R_rawTb  — framework Taubin fit on the RAW (unclustered) hits
///          R_clsTb  — framework Taubin on the arc-walk CLUSTER centroids
///        If R_rawTb is unbiased but R_cls/PRA are low, clustering pulls centroids
///        inside the arc (raises sagitta -> lowers R). If ALL are low, it's the
///        short-arc circle-fit geometry itself. Reports medians vs truth.
/// Run: root -b -q 'diag_radius_bias.C("/mnt/f/puma_sweep/output_digi_pi200_primary.root","./data/attpcsim_pi200.root",166.8)'
double medv(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
double iqrv(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

double taubinR(const std::vector<std::pair<double,double>> &pts){
   if(pts.size()<5) return -1;
   std::vector<AtHit> hs; hs.reserve(pts.size()); std::vector<const AtHit*> hp;
   for(auto&q:pts) hs.emplace_back(0, ROOT::Math::XYZPoint(q.first,q.second,0.0),1.0);
   for(auto&h:hs) hp.push_back(&h);
   AtPatterns::AtPatternCircle2D c; c.AtPattern::FitPattern(hp,-1.0);
   double R=c.GetRadius(); return (R>0&&R<1e5)?R:-1;
}

void diag_radius_bias(TString digiFile="/mnt/f/puma_sweep/output_digi_pi200_primary.root",
                      TString simFile="./data/attpcsim_pi200.root", double Rtruth=166.8, double B=4.0){
   gSystem->Load("libAtReconstruction.so");
   TFile fD(digiFile); TTree*tD=(TTree*)fD.Get("cbmsim");
   TFile fS(simFile);  TTree*tS=(TTree*)fS.Get("cbmsim");
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent",&pat);
   TClonesArray*mcP=new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint",&mcP);

   std::vector<double> rPRA,rRaw,rCls, nHitsV, arcDegV;
   for(Long64_t e=0;e<tD->GetEntries();++e){ tD->GetEntry(e); tS->GetEntry(e);
      if(!pat->GetEntries())continue;
      int nMC=mcP->GetEntries(); std::vector<double>mx(nMC),my(nMC); std::vector<int>mid(nMC);
      for(int k=0;k<nMC;++k){auto*p=(AtMCPoint*)mcP->At(k); mx[k]=p->GetX()*10; my[k]=p->GetY()*10; mid[k]=p->GetTrackID();}
      for(auto&tr:((AtPatternEvent*)pat->At(0))->GetTrackCand()){
         // truth-match to a pion (tid 0/1)
         std::map<int,int> votes; std::vector<std::pair<double,double>> raw;
         for(auto&h:tr.GetHitArray()){auto&p3=h->GetPosition(); raw.emplace_back(p3.X(),p3.Y());
            double best=9; int bid=-99; for(int k=0;k<nMC;++k){double d2=(p3.X()-mx[k])*(p3.X()-mx[k])+(p3.Y()-my[k])*(p3.Y()-my[k]); if(d2<best){best=d2;bid=mid[k];}}
            if(bid!=-99)votes[bid]++; }
         int bv=0,bt=-99; for(auto&kv:votes)if(kv.second>bv){bv=kv.second;bt=kv.first;}
         if(!(bt==0||bt==1)) continue;
         double Rp=tr.GetGeoRadius(); if(!(Rp>0&&Rp<1e5)) continue;
         // cluster centroids
         std::vector<std::pair<double,double>> cls;
         for(auto&c:*tr.GetHitClusterArray()){auto p=c.GetPosition(); cls.emplace_back(p.X(),p.Y());}
         double Rraw=taubinR(raw), Rcls=taubinR(cls);
         rPRA.push_back(Rp);
         if(Rraw>0)rRaw.push_back(Rraw);
         if(Rcls>0)rCls.push_back(Rcls);
         nHitsV.push_back(tr.GetHitArray().size());
         // arc span in degrees about the fitted centre
         auto cen=tr.GetGeoCenter(); double amin=1e9,amax=-1e9;
         for(auto&pr:raw){double a=std::atan2(pr.second-cen.second,pr.first-cen.first); amin=std::min(amin,a);amax=std::max(amax,a);}
         arcDegV.push_back((amax-amin)*180/M_PI);
      }
   }
   auto pct=[&](double R){return 100*(R-Rtruth)/Rtruth;};
   printf("\n===== radius-bias localization: %s (truth R=%.1f mm) =====\n", gSystem->BaseName(digiFile.Data()), Rtruth);
   printf("  R_PRA   (cluster-based, pipeline): median=%.1f mm (%+.1f%%)  IQR=%.1f%%  n=%zu\n",medv(rPRA),pct(medv(rPRA)),100*iqrv(rPRA)/Rtruth,rPRA.size());
   printf("  R_rawTb (Taubin on raw hits)     : median=%.1f mm (%+.1f%%)  IQR=%.1f%%  n=%zu\n",medv(rRaw),pct(medv(rRaw)),100*iqrv(rRaw)/Rtruth,rRaw.size());
   printf("  R_clsTb (Taubin on clusters)     : median=%.1f mm (%+.1f%%)  IQR=%.1f%%  n=%zu\n",medv(rCls),pct(medv(rCls)),100*iqrv(rCls)/Rtruth,rCls.size());
   printf("  median hits/track=%.0f  median arc span=%.0f deg\n", medv(nHitsV), medv(arcDegV));
   printf("===========================================================\n\n");
}
