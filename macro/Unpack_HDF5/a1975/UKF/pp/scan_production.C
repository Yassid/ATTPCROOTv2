/// @file scan_production.C
/// @brief Did every run in a production actually finish? Compare nfit against nreco, run by run.
///
/// A genfitter production can die partway through a run -- reading a ~1 GB `_reco.root` under
/// memory pressure raises `SysError ... Cannot allocate memory` -- and still leave a large,
/// perfectly readable output file behind. `refit_pp.sh` tested it with `-s`, which such a file
/// passes, so the run was marked COMPLETED. Two runs of 84 were wrong this way for eight days:
/// run_0158 (0 of 24214 events) and run_0157 (7748 of 25585).
///
/// The only reliable test is the entry count against the reco file the fit was made from. The
/// FRIB count is printed too but is NOT the reference: it legitimately runs one long on six runs.
///
/// The reco suffix is a parameter because the channels do not agree on one: the H2 chain writes
/// <run>_reco.root and the D2 chain <run>_multifit_reco.root.
///
///   root -b -q 'pp/scan_production.C'   // (p,p), the default
///   root -b -q 'pp/scan_production.C("/mnt/f/a1975/reco_pd_catima_bx/","_genfitter_pdcatbx")'
///   root -b -q 'pp/scan_production.C("/mnt/f/a1975/gf_dt_cateloss/","_multifit_genfitter_t",
///                                    "/mnt/f/a1975/reco_d2_dv1104/","_multifit_reco")'
void scan_production(TString gfDir = "/mnt/f/a1975/reco/", TString suffix = "_genfitter_pphand",
                     TString recoDir = "/mnt/f/a1975/reco/", TString recoSuffix = "_reco")
{
   auto nent = [](TString p) -> Long64_t {
      if (gSystem->AccessPathName(p)) return -2; // no file
      TFile *f = TFile::Open(p, "READ");
      if (!f || f->IsZombie()) { if (f) delete f; return -3; } // unreadable
      auto *t = (TTree *)f->Get("cbmsim");
      Long64_t n = t ? t->GetEntries() : -1;
      f->Close(); delete f;
      return n;
   };

   std::vector<TString> runs;
   void *d = gSystem->OpenDirectory(gfDir);
   const char *e;
   while ((e = gSystem->GetDirEntry(d))) {
      TString s(e);
      if (s.EndsWith(suffix + ".root")) runs.push_back(s(0, s.Index(suffix + ".root")));
   }
   gSystem->FreeDirectory(d);
   std::sort(runs.begin(), runs.end());
   if (runs.empty()) { printf("\033[1;31mno %s*%s.root found\033[0m\n", gfDir.Data(), suffix.Data()); return; }

   printf("\n%-10s %10s %10s %10s   %s\n", "run", "nfit", "nFRIB", "nreco", "verdict");
   int nBad = 0;
   Long64_t lost = 0, total = 0;
   for (auto &r : runs) {
      Long64_t nf = nent(gfDir + r + suffix + ".root");
      Long64_t nb = nent(recoDir + r + "_FRIB.root");
      Long64_t nr = nent(recoDir + r + recoSuffix + ".root");
      if (nr > 0) total += nr;
      if (nr <= 0) continue; // no reference to compare against
      if (nf == 0) {
         printf("\033[1;31m%-10s %10lld %10lld %10lld   EMPTY -- rerun this run\033[0m\n", r.Data(), nf, nb, nr);
         ++nBad; lost += nr;
      } else if (nf > 0 && nf < nr) {
         printf("\033[1;31m%-10s %10lld %10lld %10lld   TRUNCATED -- %.1f%% of the run missing, rerun it\033[0m\n",
                r.Data(), nf, nb, nr, 100.0 * (nr - nf) / nr);
         ++nBad; lost += nr - nf;
      }
   }
   if (nBad == 0) printf("  (every run complete)\n");
   printf("\n%d runs, %d incomplete; %lld of %lld reco events missing (%.2f%%)\n\n", (int)runs.size(), nBad, lost,
          total, total ? 100.0 * lost / total : 0.0);
}
