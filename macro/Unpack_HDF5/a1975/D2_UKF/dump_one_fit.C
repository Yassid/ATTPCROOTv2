/// Everything the stored fit object says about ONE track, in every production.  cmp_prod_one.C
/// printed "NO FIT" whenever ndf <= 0, which is WRONG as a label: a fit object can be present,
/// carry a vertex and a kinematics, and still have ndf <= 0.  That is a COLLAPSED fit, not an
/// absent one, and it is a different failure with a different cause.
#include <cmath>
void dump_one_fit(TString run="run_0016_multifit", Long64_t entry=1808, int tid=0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   const char *dirs[4] = {"/mnt/f/a1975/gf_dp_cateloss/", "/mnt/f/a1975/gf_dp_find/",
                          "/mnt/f/a1975/gf_dp_findscan/", "/mnt/f/a1975/gf_dp_find805/"};
   const char *nm[4] = {"baseline", "FIND coarse", "FIND 25-95", "FIND 80->5"};
   printf("\n  %s entry %lld tid %d\n", run.Data(), entry, tid);
   printf("  %-12s %6s %8s %10s %8s %7s %8s %8s %9s\n",
          "production","ncl","conv","chi2","ndf","chi2ndf","smooth","KE","vtx|xy|");
   for (int d=0; d<4; ++d) {
      TFile *f = TFile::Open(TString(dirs[d]) + run + "_genfitter_p.root");
      if (!f || f->IsZombie()) { printf("  %-12s (no file)\n", nm[d]); continue; }
      auto *t=(TTree*)f->Get("cbmsim"); TClonesArray *a=nullptr;
      t->SetBranchAddress("AtTrackingEvent",&a);
      t->GetEntry(entry);
      auto *ev=(AtTrackingEvent*)a->At(0);
      int ncl=-1;
      if(ev) for(auto &tr: ev->GetTrackArray()) if(tr.GetTrackID()==tid) ncl=tr.GetHitClusterArray()->size();
      bool any=false;
      if(ev) for(auto &ft: ev->GetFittedTracks()){
         if(!ft||ft->GetTrackID()!=tid) continue;
         auto *m = ft->GetTrackMetadata().get();
         double ndf=m->GetNdf(), chi2=m->GetChi2();
         auto k=ft->GetKinematics(); auto v=ft->GetVertex();
         printf("  %-12s %6d %8d %10.3f %8.2f %7s %8zu %8.3f %9.1f\n", nm[d], ncl,
                (int)m->GetFitConverged(), chi2, ndf,
                ndf>0?Form("%.2f",chi2/ndf):"n/a", ft->GetSmoothedPositions().size(),
                k.kineticEnergy, std::hypot(v.X(),v.Y()));
         any=true;
      }
      if(!any) printf("  %-12s %6d   -- no AtFittedTrack object for this trackID --\n", nm[d], ncl);
      f->Close();
   }
}
