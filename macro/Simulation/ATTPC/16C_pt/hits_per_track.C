/// @file hits_per_track.C
/// @brief Median hits per track candidate, for DATA or SIMULATION reco files alike.
///
/// This is the quantity the gain scan is decided on. Gain sets how many electrons reach a pad,
/// hence whether the pad crosses the PSA threshold, hence how many hits a track keeps -- which is
/// what the pattern finder and AtSpyralPID see, and therefore what an acceptance measures.
/// Choosing gain by "which value gives the highest acceptance" would be monotonic and would always
/// pick the largest; matching hits/track against the DATA is the check with a right answer.
///
/// MEDIANS, NOT MEANS: the distribution has a long tail of short fragments at one end and very
/// long beam-like tracks at the other, and a mean is dragged by both.
///
/// "NON-BEAM" IS GEOMETRIC, deliberately. Using AtSpyralPID to get a polar angle would be the
/// obvious route, but including its header alongside AtPatternEvent trips a cling assertion in
/// this ROOT build -- and more importantly the Spyral estimate is itself gain-sensitive, so
/// selecting with it would fold the thing being measured into the selection. A beam track runs
/// along the axis and stays near r = 0; anything reaching beyond rMin mm from the axis is not the
/// beam. That is a property of the hits, not of any reconstruction choice.
///
///   root -b -q 'hits_per_track.C("/mnt/f/a1975/reco/run_0106_reco.root",3000)'
/// NO explicit #include of the AT-TPC headers. The dictionaries come from the loaded libraries,
/// which is what every working macro in this tree does; including the headers by hand instead
/// trips a cling assertion ("Missing lambda call operator") in this ROOT build.
static double hpt_median(std::vector<double> v)
{
   if (v.empty()) return -1.0;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}

void hits_per_track(TString file, Long64_t maxEvt = 3000, double rMin = 20.0)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   auto *f = TFile::Open(file);
   auto *t = (f && !f->IsZombie()) ? (TTree *)f->Get("cbmsim") : nullptr;
   if (!t) { printf("  cannot open %s\n", file.Data()); return; }
   if (!t->GetBranch("AtPatternEvent")) { printf("  no AtPatternEvent in %s\n", file.Data()); return; }
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);
   std::vector<double> allH, nbH;
   const Long64_t N = std::min(maxEvt > 0 ? maxEvt : t->GetEntries(), t->GetEntries());
   long nTrk = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *p = (AtPatternEvent *)pe->At(0);
      if (!p) continue;
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         ++nTrk;
         const double nh = tr.GetHitArray().size();
         if (nh <= 0) continue;
         allH.push_back(nh);
         double rmax = 0;
         for (const auto &hit : tr.GetHitArray()) {
            if (!hit) continue;
            const auto &pos = hit->GetPosition();
            rmax = std::max(rmax, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
         }
         if (rmax > rMin) nbH.push_back(nh);
      }
   }
   printf("  %-44s ev %5lld  trk %6ld | median hits/track  ALL %6.1f (n=%zu)  NON-BEAM %6.1f (n=%zu)\n",
          gSystem->BaseName(file), N, nTrk, hpt_median(allH), allH.size(), hpt_median(nbH), nbH.size());
}
