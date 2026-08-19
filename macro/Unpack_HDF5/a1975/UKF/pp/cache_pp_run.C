/// @file cache_pp_run.C
/// @brief Cache one run's genfit proton kinematics into a flat ntuple, so the
/// interactive explorer (pp/explore_pp.C) never has to re-read the slow FRIB files.
/// Writes a TTree "pk" with: ke, theta[deg], vz[mm], chi2ndf, ic, run.
/// Ex is NOT stored (the explorer recomputes it so beam energy stays adjustable).
///
///   root -b -q 'pp/cache_pp_run.C("run_0106","/tmp/ppkin_run0106.root")'
/// (driven over all runs by xargs -P4, then hadd into /tmp/pp_kin.root)

void cache_pp_run(TString run, TString outFile, TString inDir="/mnt/f/a1975/reco/",
                  TString suffix="_genfitter_pp", int icTbLo=1000, int icTbHi=1350)
{
   gSystem->Load("libAtReconstruction.so");
   TString gf=inDir+run+suffix+".root", ff=inDir+run+"_FRIB.root";
   if(gSystem->AccessPathName(gf)||gSystem->AccessPathName(ff)){ printf("skip %s\n",run.Data()); return; }
   TFile *fg=TFile::Open(gf); TTree *tg=(TTree*)fg->Get("cbmsim");
   TFile *fc=TFile::Open(ff); TTree *tc=(TTree*)fc->Get("cbmsim");
   TClonesArray *te=nullptr,*re=nullptr;
   tg->SetBranchAddress("AtTrackingEvent",&te); tc->SetBranchAddress("AtRawEvent",&re);

   int runNo=TString(run(run.Length()-4,4)).Atoi();
   TFile out(outFile,"RECREATE");
   float ke,theta,vz,chi2ndf,ic; int rn;
   TTree *pk=new TTree("pk","genfit proton kinematics");
   pk->Branch("ke",&ke); pk->Branch("theta",&theta); pk->Branch("vz",&vz);
   pk->Branch("chi2ndf",&chi2ndf); pk->Branch("ic",&ic); pk->Branch("run",&rn);
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

   Long64_t N=std::min(tg->GetEntries(),tc->GetEntries());
   long n=0;
   for(Long64_t i=0;i<N;++i){
      tc->GetEntry(i); ic=-1;
      if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
         if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
            double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } }
      tg->GetEntry(i); if(te->GetEntries()==0) continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      for(auto&ft:ev->GetFittedTracks()){ // fitter already proton-gated + theta-windowed
         if(!ft) continue; auto&k=ft->GetKinematics(); auto&m=ft->GetTrackMetadata(); if(!m) continue;
         double ndf=m->GetNdf();
         ke=k.kineticEnergy; theta=k.theta*TMath::RadToDeg(); vz=ft->GetVertex(0).Z();
         chi2ndf=ndf>0?m->GetChi2()/ndf:1e9;
         pk->Fill(); ++n;
      }
   }
   out.cd(); pk->Write(); out.Close(); fg->Close(); fc->Close();
   printf("%s: cached %ld protons -> %s\n",run.Data(),n,outFile.Data());
}
