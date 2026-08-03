/// @file mat_audit.C
/// @brief Did the material-effects fits actually keep their material, or did they all fall back?
///
/// AtGenfitter::SetMatEffectsFallback defaults ON: a track whose material fit throws is retried
/// with material effects OFF so it still produces a result. That is sensible for yield and
/// invisible in the track count, which is exactly the problem -- running with matEffects = kTRUE
/// and reading the fitted-track count tells you nothing about whether any material was used.
/// Symptom that prompted this: matEffects on and off gave 56 vs 58 tracks in the same 3000
/// events and the SAME runtime, which material transport should not allow.
///
/// AtFitTrackMetadata records the provenance per track: GetMatEffects() is what the kept fit
/// actually used, GetMatEffectsFallback() flags that it was downgraded. Filter on these before
/// quoting any resolution from a material-effects production.
///
///   root -b -q 'mat_audit.C("/tmp/gf_kTRUE/run_0016_genfitter_t.root")'

void mat_audit(TString file)
{
   gSystem->Load("libAtReconstruction.so");
   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", file.Data()); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("no cbmsim\n"); return; }
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);

   long tot = 0, withMat = 0, fellBack = 0, good = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (te == nullptr || te->GetEntries() == 0)
         continue;
      auto *ev = dynamic_cast<AtTrackingEvent *>(te->At(0));
      if (ev == nullptr)
         continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         const auto &m = ft->GetTrackMetadata();
         if (!m)
            continue;
         ++tot;
         if (m->GetMatEffects())
            ++withMat;
         if (m->GetMatEffectsFallback())
            ++fellBack;
         if (m->GetGoodFit())
            ++good;
      }
   }
   printf("\n%s\n", gSystem->BaseName(file));
   printf("  fitted tracks          : %ld  (good fit %ld)\n", tot, good);
   printf("  kept WITH material     : %ld\n", withMat);
   printf("  downgraded to no-mat   : %ld\n", fellBack);
   if (tot > 0 && withMat == 0)
      printf("  -> EVERY track fell back: this production has no material effects in it at all\n");
}
