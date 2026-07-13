/// @file momentum_biasfix.C
/// @brief Prototype + validate the vertex-momentum bias fix. The fit reports the
///        spiral-AVERAGE momentum (~ momentum at the track mid-arc) because the
///        pion loses energy along the visible track. Recover the VERTEX momentum
///        by adding back the energy lost from the vertex to the arc mid-point
///        (~half the visible arc length) using the CATIMA dE/dx model.
/// Run: root -b -q 'momentum_biasfix.C("data/scan_p0.375.root",0.3749)'
double medOf(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
double iqrS(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

void momentum_biasfix(TString fitFile, double p0GeV, TString fitter="UKF")
{
   gSystem->Load("libAtReconstruction.so");
   const double m=139.57039, u=931.49410372;
   const double p0=p0GeV*1000.0; // MeV

   // dE/dx model identical to the one the UKF uses (CATIMA, P10 @ ~1 bar)
   auto eloss=std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
   eloss->SetProjectile(1,1,m/u);
   std::vector<std::tuple<int,int,int>> mat={{18,40,9},{6,12,1},{1,1,4}};
   eloss->SetMaterial(mat);

   TFile f(fitFile); auto*t=(TTree*)f.Get("cbmsim");
   auto*teArr=new TClonesArray("AtTrackingEvent"); t->SetBranchAddress(Form("AtTrackingEvent%s",fitter.Data()),&teArr);
   auto*patArr=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&patArr);

   std::vector<double> dpRaw, dpCorr;
   for(Long64_t e=0;e<t->GetEntries();e++){ t->GetEntry(e);
      if(teArr->GetEntries()==0||patArr->GetEntries()==0)continue;
      auto*te=(AtTrackingEvent*)teArr->At(0); auto*pe=(AtPatternEvent*)patArr->At(0);
      auto &tc=pe->GetTrackCand();
      for(const auto&ft: te->GetFittedTracks()){
         const auto&kin=ft->GetKinematics(0); double KE=kin.kineticEnergy; if(!(KE>0))continue;
         double p=std::sqrt(KE*KE+2*KE*m);           // measured (mid-arc) momentum MeV
         // visible arc length from the matching PRA track: R * angular span
         int tid=ft->GetTrackID(); if(tid<0||tid>=(int)tc.size())continue;
         auto &tr=tc[tid]; double R=tr.GetGeoRadius(); auto ctr=tr.GetGeoCenter();
         double phiMin=1e9,phiMax=-1e9;
         for(const auto&h: tr.GetHitArray()){ const auto&pp=h->GetPosition();
            double ph=std::atan2(pp.Y()-ctr.second, pp.X()-ctr.first);
            phiMin=std::min(phiMin,ph); phiMax=std::max(phiMax,ph); }
         double dPhi=phiMax-phiMin; if(dPhi>M_PI)dPhi=2*M_PI-dPhi;
         double arc=(R>0&&R<1e4)?R*dPhi:0;           // visible arc length [mm]
         // energy lost from vertex to arc mid-point ~ half the visible arc
         double halfArc=0.5*arc;
         double KEc=KE, pc=p;
         if(halfArc>0){ double dedx=eloss->GetdEdx(KE); double eLost=dedx*halfArc;
            KEc=KE+eLost; pc=std::sqrt(KEc*KEc+2*KEc*m); }
         dpRaw.push_back(100*(p-p0)/p0);
         dpCorr.push_back(100*(pc-p0)/p0);
      }
   }
   printf("  p0=%.0f MeV  n=%zu :  RAW bias %+.1f%% (sig %.1f%%)  ->  CORRECTED bias %+.1f%% (sig %.1f%%)\n",
          p0, dpRaw.size(), medOf(dpRaw), iqrS(dpRaw), medOf(dpCorr), iqrS(dpCorr));
}
