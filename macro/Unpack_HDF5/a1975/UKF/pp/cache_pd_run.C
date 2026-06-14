/// @file cache_pd_run.C
/// @brief Cache one run's genfit DEUTERON kinematics into a flat ntuple for the
/// interactive 16C(p,d)15C explorer (pp/explore_pd.C). Sibling of cache_pp_run.C.
/// Writes a TTree "pk" with: ke, theta[deg], vz[mm], chi2ndf, ic, run. Ex is NOT
/// stored (the explorer recomputes it so beam energy stays adjustable).
///
/// The deuteron fits live in gfDir (reco_pd/, *_genfitter_pd.root) while the FRIB
/// ion-chamber files live in fribDir (reco/, *_FRIB.root).
///
///   root -b -q 'pp/cache_pd_run.C("run_0106","/tmp/pdkin_run0106.root")'
/// (driven over all runs by xargs -P4, then hadd into /tmp/pd_kin.root)

void cache_pd_run(TString run, TString outFile, TString gfDir="/mnt/f/a1975/reco_pd/",
                  TString fribDir="/mnt/f/a1975/reco/", TString suffix="_genfitter_pd",
                  int icTbLo=1000, int icTbHi=1350)
{
   gSystem->Load("libAtReconstruction.so");
   TString gf=gfDir+run+suffix+".root", ff=fribDir+run+"_FRIB.root";
   if(gSystem->AccessPathName(gf)){ printf("skip %s (no %s)\n",run.Data(),gf.Data()); return; }
   bool haveFrib = !gSystem->AccessPathName(ff);
   TFile *fg=TFile::Open(gf); TTree *tg=(TTree*)fg->Get("cbmsim");
   TFile *fc=haveFrib?TFile::Open(ff):nullptr; TTree *tc=haveFrib?(TTree*)fc->Get("cbmsim"):nullptr;
   TClonesArray *te=nullptr,*re=nullptr;
   tg->SetBranchAddress("AtTrackingEvent",&te);
   if(tc) tc->SetBranchAddress("AtRawEvent",&re);

   int runNo=TString(run(run.Length()-4,4)).Atoi();
   TFile out(outFile,"RECREATE");
   float ke,theta,vz,chi2ndf,ic; int rn;
   TTree *pk=new TTree("pk","genfit deuteron kinematics");
   pk->Branch("ke",&ke); pk->Branch("theta",&theta); pk->Branch("vz",&vz);
   pk->Branch("chi2ndf",&chi2ndf); pk->Branch("ic",&ic); pk->Branch("run",&rn);
   rn=runNo;

   Long64_t N = tc ? std::min(tg->GetEntries(),tc->GetEntries()) : tg->GetEntries();
   long n=0;
   for(Long64_t i=0;i<N;++i){
      ic=-1;
      if(tc){ tc->GetEntry(i);
         if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
            if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
               double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } } }
      tg->GetEntry(i); if(te->GetEntries()==0) continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      for(auto&ft:ev->GetFittedTracks()){ // fitter already deuteron-gated + theta-windowed
         if(!ft) continue; auto&k=ft->GetKinematics(); auto&m=ft->GetTrackMetadata(); if(!m) continue;
         double ndf=m->GetNdf();
         ke=k.kineticEnergy; theta=k.theta*TMath::RadToDeg(); vz=ft->GetVertex(0).Z();
         chi2ndf=ndf>0?m->GetChi2()/ndf:1e9;
         pk->Fill(); ++n;
      }
   }
   out.cd(); pk->Write(); out.Close(); fg->Close(); if(fc) fc->Close();
   printf("%s: cached %ld deuterons -> %s%s\n",run.Data(),n,outFile.Data(), haveFrib?"":" (no IC)");
}
