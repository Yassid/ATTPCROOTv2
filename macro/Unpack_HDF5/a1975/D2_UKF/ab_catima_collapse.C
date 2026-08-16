/// @file ab_catima_collapse.C
/// @brief Collapse rate (ndf <= 0) vs KE_Brho, per arm.
///
/// KE_Brho is the Spyral arc fit and is independent of genfit, so it classifies a track even
/// when the genfit fit collapsed -- which is the only way to ask WHICH tracks are being lost.
///
/// The question this answers: the dE/dx table serves beta*gamma < 0.05 (KE < 3.5 MeV for a
/// triton), yet removing it collapses ~48% of tracks while only ~18% of tracks START below
/// 3.5 MeV. If the loss is spread across all energies rather than confined to the low bin,
/// that is because a track SLOWS THROUGH the dead region on its way to stopping, so the table
/// matters for most of the sample and not just for the low branch.

#include <cmath>
#include <cstdio>
#include <vector>

void ab_catima_collapse(TString dir = "/mnt/f/a1975/")
{
   gSystem->Load("libAtReconstruction.so");
   const char *arms[] = {"nomat", "off", "on", "nostrag", "notable"};
   const double edges[] = {0, 2, 3.5, 6, 10, 20, 1e9};
   const int nb = 6;

   printf("\ncollapse rate (ndf<=0 or chi2<=0) by KE_Brho at the FIRST measurement point\n");
   printf("%-9s", "arm");
   for (int b = 0; b < nb; ++b)
      printf("  %10s", Form("%g-%g", edges[b], edges[b + 1] > 1e8 ? 999.0 : edges[b + 1]));
   printf("  %10s\n", "all");

   for (auto arm : arms) {
      TString f = dir + "gf_dt_ab_" + arm + "/run_0031_multifit_genfitter_t.root";
      TFile *fl = TFile::Open(f);
      if (!fl || fl->IsZombie()) {
         printf("%-9s  (missing)\n", arm);
         continue;
      }
      TTree *t = (TTree *)fl->Get("cbmsim");
      TClonesArray *te = nullptr, *pe = nullptr;
      t->SetBranchAddress("AtTrackingEvent", &te);
      t->SetBranchAddress("AtPIDEvent", &pe);
      std::vector<long> n(nb, 0), bad(nb, 0);
      long nAll = 0, badAll = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!te || !pe || te->GetEntries() == 0 || pe->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         auto *pidev = (AtPIDEvent *)pe->At(0);
         if (!ev || !pidev)
            continue;
         std::map<int, AtFittedTrack *> fmap;
         for (auto &ft : ev->GetFittedTracks())
            if (ft)
               fmap[ft->GetTrackID()] = ft.get();
         for (auto &sr : pidev->GetSpyral()) {
            if (!sr.valid || sr.brho <= 0)
               continue;
            auto it = fmap.find(sr.trackID);
            if (it == fmap.end())
               continue;
            const auto &m = it->second->GetTrackMetadata();
            if (!m)
               continue;
            const double p = 299.792458 * sr.brho, mt = 2808.921;
            const double ke = std::sqrt(p * p + mt * mt) - mt;
            const bool isBad = !(m->GetNdf() > 0 && m->GetChi2() > 0);
            for (int b = 0; b < nb; ++b)
               if (ke >= edges[b] && ke < edges[b + 1]) {
                  n[b]++;
                  bad[b] += isBad;
               }
            nAll++;
            badAll += isBad;
         }
      }
      printf("%-9s", arm);
      for (int b = 0; b < nb; ++b)
         printf("  %9s", n[b] ? Form("%.0f%%(%ld)", 100.0 * bad[b] / n[b], n[b]) : "-");
      printf("  %9s\n", nAll ? Form("%.0f%%(%ld)", 100.0 * badAll / nAll, nAll) : "-");
      fl->Close();
   }
   printf("\n");
}
