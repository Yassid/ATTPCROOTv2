/// How many stored fits are COLLAPSED -- a fit object exists but ndf <= 0, so chi2/ndf does not
/// exist and every downstream cut on it is silently skipped.  Counted per production, because the
/// question is whether truncation created them.
#include <cmath>
void collapse_count_dp(TString gfDir, TString runs="run_0016_multifit,run_0020_multifit")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TObjArray *ra = runs.Tokenize(",");
   long nFit=0, nColl=0, nNotConv=0, nOnePt=0;
   std::vector<double> prevSig;
   for (int i=0;i<ra->GetEntries();++i){
      TString rt=((TObjString*)ra->At(i))->GetString();
      TFile *f=TFile::Open(gfDir+rt+"_genfitter_p.root");
      if(!f||f->IsZombie()) continue;
      auto *t=(TTree*)f->Get("cbmsim"); TClonesArray*a=nullptr;
      t->SetBranchAddress("AtTrackingEvent",&a);
      for(Long64_t e=0;e<t->GetEntries();++e){
         t->GetEntry(e);
         if(a->GetEntries()==0) continue;
         auto *ev=(AtTrackingEvent*)a->At(0); if(!ev) continue;
         std::vector<double> sig;
         for(auto &tr: ev->GetTrackArray()) sig.push_back(tr.GetHitClusterArray()->size());
         if(!sig.empty()&&sig==prevSig) continue;      // duplicate entry guard
         prevSig=sig;
         for(auto &ft: ev->GetFittedTracks()){
            if(!ft) continue;
            ++nFit;
            auto *m=ft->GetTrackMetadata().get();
            if(!(m->GetNdf()>0)) ++nColl;
            if(!m->GetFitConverged()) ++nNotConv;
            if(ft->GetSmoothedPositions().size()<=1) ++nOnePt;
         }
      }
      f->Close();
   }
   printf("  %-34s fits %6ld | ndf<=0 %5ld (%4.1f%%) | not converged %6ld (%4.1f%%) | <=1 stored pos %5ld (%4.1f%%)\n",
          gSystem->BaseName(gfDir.Data()), nFit, nColl, nFit?100.0*nColl/nFit:0.0,
          nNotConv, nFit?100.0*nNotConv/nFit:0.0, nOnePt, nFit?100.0*nOnePt/nFit:0.0);
}
