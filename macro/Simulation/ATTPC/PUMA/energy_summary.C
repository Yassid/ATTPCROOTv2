/// @file energy_summary.C
/// @brief UKF vs genfit efficiency + resolution across pion momenta (150-500 MeV/c),
///        on the production pipeline (primaryOnly + CircleMerge). Per (energy, fitter):
///          efficiency  = distinct reconstructable pions with a valid matched fit / all
///                        reconstructable pions (>= NMIN gas MC points)
///          resolution  = IQR of Delta p/p (vs vertex momentum) for matched fits
///          bias        = median Delta p/p  (includes the Cu/Al material offset)
///          dtheta      = IQR of theta - 90 deg
///        Writes a table and two plots (resolution & efficiency vs momentum).
/// Run: root -b -q energy_summary.C
double medv(std::vector<double>v){if(v.empty())return 0;std::sort(v.begin(),v.end());return v[v.size()/2];}
double iqrv(std::vector<double>v){if(v.size()<4)return 0;std::sort(v.begin(),v.end());return(v[3*v.size()/4]-v[v.size()/4])/1.349;}

struct Pt { const char* tag; double p0; };

struct Res { double eff, sig, med, dth; int nfit, nreco; };

Res analyze(TString digi, TString sim, double p0, int fitter /*0=UKF,1=genfit*/, int NMIN=15, double B=4.0){
   const double m=139.57039;
   Res r{0,0,0,0,0,0};
   TFile fD(digi); TTree*tD=(TTree*)fD.Get("cbmsim");
   TFile fS(sim);  TTree*tS=(TTree*)fS.Get("cbmsim");
   if(!tD||!tS) return r;
   const char* br = fitter==0?"AtTrackingEventUKF":"AtTrackingEventGenfit";
   TClonesArray*fe=new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress(br,&fe);
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent",&pat);
   TClonesArray*mcP=new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint",&mcP);
   std::vector<double> dp, dth; int nreco=0; std::set<std::pair<Long64_t,int>> fittedPions;
   Long64_t nE=std::min(tD->GetEntries(),tS->GetEntries());
   for(Long64_t e=0;e<nE;++e){ tD->GetEntry(e); tS->GetEntry(e);
      int nMC=mcP->GetEntries(); std::vector<double>mx(nMC),my(nMC); std::vector<int>mid(nMC);
      int cnt[2]={0,0};
      for(int k=0;k<nMC;++k){auto*p=(AtMCPoint*)mcP->At(k); mx[k]=p->GetX()*10; my[k]=p->GetY()*10; mid[k]=p->GetTrackID();
         if(mid[k]==0||mid[k]==1)cnt[mid[k]]++; }
      for(int id=0;id<2;++id) if(cnt[id]>=NMIN) nreco++; // reconstructable pions this event
      if(!pat->GetEntries()||!fe->GetEntries()) continue;
      auto&tc=((AtPatternEvent*)pat->At(0))->GetTrackCand();
      for(const auto&ft:((AtTrackingEvent*)fe->At(0))->GetFittedTracks()){
         double KE=ft->GetKinematics(0).kineticEnergy; if(!(KE>0))continue;
         int tid=ft->GetTrackID(); if(tid<0||tid>=(int)tc.size())continue;
         // truth-match the PRA track's hits to a pion
         std::map<int,int> votes;
         for(auto&h:tc[tid].GetHitArray()){auto&p3=h->GetPosition(); double best=9;int bid=-99;
            for(int k=0;k<nMC;++k){double d2=(p3.X()-mx[k])*(p3.X()-mx[k])+(p3.Y()-my[k])*(p3.Y()-my[k]); if(d2<best){best=d2;bid=mid[k];}}
            if(bid!=-99)votes[bid]++; }
         int bv=0,bt=-99; for(auto&kv:votes)if(kv.second>bv){bv=kv.second;bt=kv.first;}
         if(!(bt==0||bt==1)) continue;
         double p=sqrt(KE*KE+2*KE*m);
         dp.push_back(100*(p-p0)/p0);
         dth.push_back(ft->GetKinematics(0).theta*180/M_PI-90.0);
         fittedPions.insert({e,bt});
      }
   }
   r.nreco=nreco; r.nfit=dp.size();
   r.eff = nreco? 100.0*fittedPions.size()/nreco : 0;
   r.sig = iqrv(dp); r.med = medv(dp); r.dth = iqrv(dth);
   return r;
}

void energy_summary(){
   gSystem->Load("libAtReconstruction.so"); gStyle->SetOptStat(0);
   std::vector<Pt> pts = {{"p150",150},{"p200",200},{"p300",300},{"p375",375},{"p500",500}};
   TString D="/mnt/f/puma_sweep/", S="./data/";
   printf("\n============ UKF vs genfit — efficiency + resolution vs pion momentum ============\n");
   printf(" primaryOnly + CircleMerge; eff = pions fitted / reconstructable (>=15 gas pts)\n");
   printf(" %-6s | %-30s | %-30s\n","p MeV","UKF   eff  sigp  bias  dth","genfit eff  sigp  bias  dth");
   printf(" -------------------------------------------------------------------------------------\n");
   std::vector<double> vp, effU,sigU,effG,sigG, medU,medG;
   for(auto&P:pts){
      TString dg=D+"output_digi_es_"+P.tag+".root", sm=S+"attpcsim_"+P.tag+".root";
      Res u=analyze(dg,sm,P.p0,0), g=analyze(dg,sm,P.p0,1);
      printf(" %-6.0f | %5.0f%% %5.1f%% %+5.1f%% %5.2f | %5.0f%% %5.1f%% %+5.1f%% %5.2f\n",
             P.p0, u.eff,u.sig,u.med,u.dth,  g.eff,g.sig,g.med,g.dth);
      vp.push_back(P.p0); effU.push_back(u.eff); sigU.push_back(u.sig); medU.push_back(u.med);
      effG.push_back(g.eff); sigG.push_back(g.sig); medG.push_back(g.med);
   }
   printf(" -------------------------------------------------------------------------------------\n");
   printf(" sigp=IQR resolution; bias=median (includes Cu/Al material offset, grows at low p)\n\n");

   auto mk=[&](std::vector<double>&y,int col,int mst){auto*gr=new TGraph(vp.size(),vp.data(),y.data());
      gr->SetLineColor(col);gr->SetMarkerColor(col);gr->SetMarkerStyle(mst);gr->SetLineWidth(2);gr->SetMarkerSize(1.5);return gr;};
   auto*c=new TCanvas("es","energy scan",1300,560); c->Divide(2,1);
   c->cd(1); gPad->SetGrid();
   auto*gu=mk(sigU,kBlue+1,20),*gg=mk(sigG,kRed+1,21);
   auto*mg=new TMultiGraph(); mg->Add(gu);mg->Add(gg);
   mg->SetTitle("momentum resolution vs p;pion momentum [MeV/c];#sigma_{p}/p (IQR) [%]");
   mg->Draw("ALP"); mg->SetMinimum(0);
   auto*l1=new TLegend(0.55,0.72,0.88,0.88); l1->AddEntry(gu,"UKF","lp"); l1->AddEntry(gg,"genfit","lp"); l1->Draw();
   c->cd(2); gPad->SetGrid();
   auto*eu=mk(effU,kBlue+1,20),*eg=mk(effG,kRed+1,21);
   auto*mg2=new TMultiGraph(); mg2->Add(eu);mg2->Add(eg);
   mg2->SetTitle("efficiency vs p;pion momentum [MeV/c];efficiency [%]");
   mg2->Draw("ALP"); mg2->SetMinimum(0); mg2->SetMaximum(105);
   auto*l2=new TLegend(0.55,0.20,0.88,0.36); l2->AddEntry(eu,"UKF","lp"); l2->AddEntry(eg,"genfit","lp"); l2->Draw();
   c->SaveAs("./data/energy_summary.png");
   printf(" wrote ./data/energy_summary.png\n\n");
}
