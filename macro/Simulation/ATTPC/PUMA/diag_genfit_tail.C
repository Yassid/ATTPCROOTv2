/// @file diag_genfit_tail.C
/// @brief Characterise the genfit vertex dr-tail on the branch-8 sample: split
///        tracks into core (dr<20 mm) vs tail (dr>=20 mm) and compare fit quality,
///        charge correctness, extrapolation length, and momentum.
/// Run: root -b -q diag_genfit_tail.C
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void diag_genfit_tail(TString digiFile="./data/output_digi_both8.root", TString simFile="./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   const double m_pi=139.57039, p0=374.9, kTol=3.0, drCut=20.0;
   TFile fD(digiFile); TTree*tD=(TTree*)fD.Get("cbmsim");
   TFile fS(simFile);  TTree*tS=(TTree*)fS.Get("cbmsim");
   TClonesArray*gf=new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress("AtTrackingEventGenfit",&gf);
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent",&pat);
   TClonesArray*mcPts=new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint",&mcPts);
   TClonesArray*mcTrk=new TClonesArray("AtMCTrack"); tS->SetBranchAddress("MCTrack",&mcTrk);

   struct Grp{ int n=0,wrongQ=0,qTot=0,conv=0; std::vector<double> extrap,chi2ndf,dp,firstR; };
   Grp core,tail;
   Long64_t nE=std::min(tD->GetEntries(),tS->GetEntries());
   for(Long64_t e=0;e<nE;++e){ tD->GetEntry(e); tS->GetEntry(e);
      if(!pat->GetEntries()||!gf->GetEntries())continue;
      auto*patE=(AtPatternEvent*)pat->At(0); auto&tracks=patE->GetTrackCand();
      int nMC=mcPts->GetEntries(); std::vector<double>mx(nMC),my(nMC); std::vector<int>mpdg(nMC);
      for(int k=0;k<nMC;++k){auto*mp=(AtMCPoint*)mcPts->At(k); mx[k]=mp->GetX()*10; my[k]=mp->GetY()*10;
         int t=mp->GetTrackID(); auto*mt=(t>=0&&t<mcTrk->GetEntries())?(AtMCTrack*)mcTrk->At(t):nullptr; mpdg[k]=mt?mt->GetPdgCode():0;}
      auto*te=(AtTrackingEvent*)gf->At(0);
      for(auto&ft:te->GetFittedTracks()){
         auto&pr=ft->GetTrackPropertiesStruct();
         double dr=std::hypot(pr.initialPositionXtr.X(),pr.initialPositionXtr.Y());
         double firstR=std::hypot(pr.initialPosition.X(),pr.initialPosition.Y());
         Grp&g=(dr<drCut)?core:tail; g.n++;
         g.extrap.push_back(pr.extrapolatedDistance);
         g.firstR.push_back(firstR);
         double KE=ft->GetKinematics(0).kineticEnergy; double p=std::sqrt(KE*KE+2*KE*m_pi);
         g.dp.push_back((p-p0)/p0);
         auto&meta=ft->GetTrackMetadata();
         if(meta){ double c=meta->GetChi2(); int nd=meta->GetNdf(); g.chi2ndf.push_back(nd>0?c/nd:1e9); if(meta->GetFitConverged())g.conv++; }
         // charge correctness
         auto&pinfo=ft->GetParticleInfo(0); int fitSign=pinfo.charge>0?1:(pinfo.charge<0?-1:0);
         int tid=ft->GetTrackID();
         if(fitSign&&tid>=0&&tid<(int)tracks.size()){ std::map<int,int>votes;
            for(auto&h:tracks[tid].GetHitArray()){auto pp=h->GetPosition(); double best=kTol*kTol; int bp=0;
               for(int k=0;k<nMC;++k){double d2=(pp.X()-mx[k])*(pp.X()-mx[k])+(pp.Y()-my[k])*(pp.Y()-my[k]); if(d2<best){best=d2;bp=mpdg[k];}}
               if(bp)votes[bp]++;}
            int tp=0,bv=0; for(auto&kv:votes)if(kv.second>bv){bv=kv.second;tp=kv.first;}
            if(tp){g.qTot++; int ts=tp>0?1:-1; if(ts!=fitSign)g.wrongQ++;}}
      }
   }
   auto rep=[&](const char*nm,Grp&g){
      printf("\n%s : %d tracks\n",nm,g.n); if(!g.n)return;
      printf("  converged      : %.0f%%\n",100.0*g.conv/g.n);
      printf("  chi2/ndf median: %.2f\n",med(g.chi2ndf));
      printf("  extrapLen med  : %.1f mm\n",med(g.extrap));
      printf("  firstClusterR  : %.1f mm\n",med(g.firstR));
      printf("  p bias median  : %+.1f%%\n",100*med(g.dp));
      printf("  wrong charge   : %.1f%% (%d/%d)\n",g.qTot?100.0*g.wrongQ/g.qTot:0.,g.wrongQ,g.qTot);
   };
   printf("######## genfit dr-tail diagnostic (drCut=%.0f mm) ########",drCut);
   rep("CORE (dr<20mm)",core);
   rep("TAIL (dr>=20mm)",tail);
   printf("\ntail fraction: %.1f%% (%d/%d)\n",100.0*tail.n/(core.n+tail.n),tail.n,core.n+tail.n);
}
