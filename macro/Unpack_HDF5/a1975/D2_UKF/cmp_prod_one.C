/// What each of the four (d,p) productions did with ONE track.  The GUI is pointed at
/// gf_dp_cateloss, which is the BASELINE and truncates nothing, so "the adopted fit" seen there
/// has by construction never had a prefix tried on it.  This says whether the FIND productions
/// did any better, and what they kept.
#include <cmath>
void cmp_prod_one(TString run = "run_0020_multifit", Long64_t entry = 1597, int tid = 0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   const char *dirs[4] = {"/mnt/f/a1975/gf_dp_cateloss/", "/mnt/f/a1975/gf_dp_find/",
                          "/mnt/f/a1975/gf_dp_findscan/", "/mnt/f/a1975/gf_dp_find805/"};
   const char *nm[4] = {"baseline (no FIND)", "FIND 90/75/50/25", "FIND 5% 25-95", "FIND 5% 80->5"};
   printf("\n  %s entry %lld tid %d\n", run.Data(), entry, tid);
   printf("  %-20s %8s %9s %10s %8s %8s %9s\n",
          "production", "ncl", "KE", "theta", "chi2ndf", "fitpts", "vtx |xy|");
   for (int d = 0; d < 4; ++d) {
      TFile *f = TFile::Open(TString(dirs[d]) + run + "_genfitter_p.root");
      if (!f || f->IsZombie()) { printf("  %-20s  (no file)\n", nm[d]); continue; }
      auto *t = (TTree *)f->Get("cbmsim");
      TClonesArray *a = nullptr; t->SetBranchAddress("AtTrackingEvent", &a);
      if (entry >= t->GetEntries()) { printf("  %-20s  (entry out of range)\n", nm[d]); f->Close(); continue; }
      t->GetEntry(entry);
      auto *ev = (AtTrackingEvent *)a->At(0);
      int ncl = -1;
      if (ev) for (auto &tr : ev->GetTrackArray())
         if (tr.GetTrackID() == tid) ncl = tr.GetHitClusterArray()->size();
      bool found = false;
      if (ev) for (auto &ft : ev->GetFittedTracks()) {
         if (!ft || ft->GetTrackID() != tid) continue;
         double ndf = ft->GetTrackMetadata()->GetNdf();
         double c2 = (ndf > 0) ? ft->GetTrackMetadata()->GetChi2()/ndf : -1;
         auto k = ft->GetKinematics();
         auto v = ft->GetVertex();
         // fit points is a LOWER BOUND: it counts stored smoothed positions, and a fit can store
         // fewer than it used -- see the caveat in project_a1975_dp_truncation
         printf("  %-20s %8d %9.3f %10.2f %8.2f %8zu %9.1f\n", nm[d], ncl,
                k.kineticEnergy, k.theta*TMath::RadToDeg(), c2,
                ft->GetSmoothedPositions().size(), std::hypot(v.X(), v.Y()));
         found = true;
      }
      if (!found) printf("  %-20s %8d %9s %10s %8s %8s %9s\n", nm[d], ncl, "-", "-", "NO FIT", "-", "-");
      f->Close();
   }
   printf("\n  FIND mode 0 keeps the full-length fit UNLESS chi2/ndf < c2Max (0.1) fails, and if no\n"
          "  prefix reaches 0.1 it RESTORES the full-length one -- so a truncated fit that is much\n"
          "  better but still above 0.1 is computed and then thrown away.\n");
}
