/// @file cache_pd_run.C
/// @brief Cache one run's genfit DEUTERON kinematics into a flat ntuple for the
/// interactive 16C(p,d)15C explorer (pp/explore_pd.C). Sibling of cache_pp_run.C.
/// Writes a TTree "pk" with: ke, theta[deg], vz[mm], chi2ndf, ic, run, kefit, thetafit. Ex is
/// NOT stored (the explorer recomputes it so beam energy stays adjustable). ke/theta are the
/// back-extrapolated (corrected) slot and kefit/thetafit the raw fit -- see the branch block.
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
   float ke,theta,vz,chi2ndf,ic,kefit,thetafit; int rn;
   TTree *pk=new TTree("pk","genfit deuteron kinematics");
   pk->Branch("ke",&ke); pk->Branch("theta",&theta); pk->Branch("vz",&vz);
   pk->Branch("chi2ndf",&chi2ndf); pk->Branch("ic",&ic); pk->Branch("run",&rn);
   // BOTH kinematics slots, same names and same meaning as the (d,t) caches and as mkexp_pp
   // documents them, so the explorer can select either column without deriving a second file:
   //   ke/theta       GetKinematicsXtr() -- back-extrapolated to the beam axis, the CORRECTED value
   //   kefit/thetafit GetKinematics()    -- the raw fit at the first measurement point
   // This macro previously stored GetKinematics() as `ke`. For every production built without
   // back-extrapolation the two slots are bit-identical, so no existing cache changes meaning
   // when it is rebuilt; the columns only diverge once a fit runs with backExtrap AND matEffects.
   pk->Branch("kefit",&kefit); pk->Branch("thetafit",&thetafit);
   rn=runNo;

   // ENTRY-COUNT CHECK. N below is min(fit, FRIB), so a SHORT FRIB file does not fail -- it
   // silently truncates the run and still prints the same "cached N" line as a healthy one.
   // run_0148's FRIB tree holds 1 entry against 23514 in reco, so this loop cached ONE event of
   // it: the run is absent from the caches entirely and nothing said so. It surfaced only when a
   // later macro happened to print entry counts.
   // A LONGER FRIB tree is harmless: six of the H2 runs end with a junk entry carrying event ID
   // -1, and because it sits at the END the indices still align for every real event.
   {
      const Long64_t nfit = tg->GetEntries(), nfrib = tc ? tc->GetEntries() : -1;
      if (nfrib >= 0 && nfrib < nfit)
         printf("\033[1;31m%s: FRIB SHORT -- %lld entries against %lld in the fit file; %lld events "
                "(%.1f%%) of this run are being DROPPED\033[0m\n",
                run.Data(), nfrib, nfit, nfit - nfrib, nfit ? 100.0 * (nfit - nfrib) / nfit : 0.0);
      else if (nfrib > nfit)
         printf("%s: FRIB has %lld entries against %lld -- trailing junk, harmless\n", run.Data(), nfrib, nfit);
   }

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
         if(!ft) continue; auto&k=ft->GetKinematicsXtr(); auto&kr=ft->GetKinematics();
         auto&m=ft->GetTrackMetadata(); if(!m) continue;
         double ndf=m->GetNdf();
         ke=k.kineticEnergy; theta=k.theta*TMath::RadToDeg(); vz=ft->GetVertex(0).Z();
         kefit=kr.kineticEnergy; thetafit=kr.theta*TMath::RadToDeg();
         chi2ndf=ndf>0?m->GetChi2()/ndf:1e9;
         pk->Fill(); ++n;
      }
   }
   out.cd(); pk->Write(); out.Close(); fg->Close(); if(fc) fc->Close();
   printf("%s: cached %ld deuterons -> %s%s\n",run.Data(),n,outFile.Data(), haveFrib?"":" (no IC)");
}
