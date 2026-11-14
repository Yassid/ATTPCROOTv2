/// @file dump_pd_tracks.C  (p,d) tracks of one run with their reco ENTRY index, for joining to raw FRIB mesh features.
/// Same loop as a1975_analysis/pd/production/cache_pd_run.C on the adopted reco_pd_catima_bx (fitter already deuteron-gated), plus the entry.
/// Columns: run entry ke theta vz chi2ndf ic(FRIB.root, TB 1000-1350)
void dump_pd_tracks(TString run = "run_0106", TString out = "/mnt/f/a1975/sim_rebuild/trig/pd/tracks_0106.txt",
                    TString gfDir = "/mnt/f/a1975/reco_pd_catima_bx/", TString fribDir = "/mnt/f/a1975/reco/", TString suffix = "_genfitter_pdcatbx")
{
   gErrorIgnoreLevel = kError; gSystem->Load("libAtReconstruction.so");
   TFile *fg = TFile::Open(gfDir + run + suffix + ".root"); TTree *tg = (TTree *)fg->Get("cbmsim");
   TString ff = fribDir + run + "_FRIB.root"; TFile *fc = gSystem->AccessPathName(ff) ? nullptr : TFile::Open(ff); TTree *tc = fc ? (TTree *)fc->Get("cbmsim") : nullptr;
   TClonesArray *te = nullptr, *re = nullptr; tg->SetBranchAddress("AtTrackingEvent", &te); if (tc) tc->SetBranchAddress("AtRawEvent", &re);
   int runNo = TString(run(run.Length() - 4, 4)).Atoi(); FILE *o = fopen(out, "w"); long n = 0;
   Long64_t N = tc ? std::min(tg->GetEntries(), tc->GetEntries()) : tg->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      double ic = -1;
      if (tc) { tc->GetEntry(i); if (re->GetEntries() > 0) { auto *raw = (AtRawEvent *)re->At(0); if (raw && !raw->GetGenTraces().empty()) { auto &adc = raw->GetGenTraces()[0]->GetADC(); double mx = -1e9; for (int b = 1000; b < 1350 && b < (int)adc.size(); ++b) mx = std::max(mx, (double)adc[b]); ic = mx; } } }
      tg->GetEntry(i); if (te->GetEntries() == 0) continue; auto *ev = (AtTrackingEvent *)te->At(0); if (!ev) continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue; auto &k = ft->GetKinematicsXtr(); auto &m = ft->GetTrackMetadata(); if (!m) continue; double ndf = m->GetNdf();
         fprintf(o, "%d %lld %.4f %.4f %.2f %.4f %.1f\n", runNo, i, k.kineticEnergy, k.theta * TMath::RadToDeg(), ft->GetVertex(0).Z(), ndf > 0 ? m->GetChi2() / ndf : 1e9, ic); ++n;
      }
   }
   fclose(o); printf("done %s %ld\n", run.Data(), n);
}
