/// @file compare_reco_3Hed.C
/// @brief Hit-by-hit diff of two reconstruction outputs. Used to prove that a code change left
///        normal-mode running untouched.
///
///   root -b -q 'compare_reco_3Hed.C("reference_reco.root","new_reco.root")'
///
/// THIS IS ONLY A REGRESSION TEST BECAUSE THE DIGITISATION IS REPRODUCIBLE. Nothing in the
/// simulation-reconstruction chain calls SetSeed, so gRandom keeps ROOT's default TRandom3
/// sequence and two runs over the same input draw identical numbers -- including the diffusion
/// throws and the gain fluctuations. If a seed is ever introduced, or if the number of gRandom
/// calls per event changes, this comparison stops being a regression test and starts being a
/// statistical one, and it will report spurious differences. It does not try to detect that; the
/// warning is here.
///
/// Compares the AtEvent hit collections: per-event hit COUNT first, then every hit's position and
/// charge. Reports the first disagreement with enough context to find it, rather than only a
/// verdict.
void compare_reco_3Hed(TString refFile, TString newFile, Double_t tol = 1e-6, Int_t maxReport = 5)
{
   gSystem->Load("libAtReconstruction.so");

   TFile *fA = TFile::Open(refFile), *fB = TFile::Open(newFile);
   if (!fA || fA->IsZombie() || !fB || fB->IsZombie()) {
      printf("cannot open both files\n");
      return;
   }
   TTree *tA = (TTree *)fA->Get("cbmsim"), *tB = (TTree *)fB->Get("cbmsim");
   // REFUSE A TRUNCATED FILE INSTEAD OF CRASHING ON IT. A killed reconstruction leaves a large,
   // non-empty file whose TTree was never written; ROOT "recovers" a few keys, cbmsim is not among
   // them, and dereferencing the null tree segfaults -- which is what happened the first time this
   // ran. An existence-or-size check upstream cannot see this; only reading the product can.
   if (!tA || !tB) {
      printf("\n  [DIFFER] cannot read a cbmsim tree from %s -- the file is truncated or was never\n"
             "  closed (a killed job leaves exactly this). Nothing compared.\n\n",
             !tA ? refFile.Data() : newFile.Data());
      return;
   }
   TClonesArray *eA = nullptr, *eB = nullptr;
   tA->SetBranchAddress("AtEventH", &eA);
   tB->SetBranchAddress("AtEventH", &eB);

   const Long64_t NA = tA->GetEntries(), NB = tB->GetEntries();
   const Long64_t N = std::min(NA, NB);
   printf("\n=== reco comparison ===\n  reference : %s (%lld events)\n  new       : %s (%lld events)\n",
          refFile.Data(), NA, newFile.Data(), NB);
   if (NA != NB)
      printf("  *** DIFFERENT EVENT COUNTS *** comparing the first %lld\n", N);

   long nHitA = 0, nHitB = 0, nCountMismatch = 0, nPosMismatch = 0, reported = 0;
   double worst = 0;

   for (Long64_t i = 0; i < N; ++i) {
      tA->GetEntry(i);
      tB->GetEntry(i);
      auto *evA = eA->GetEntriesFast() ? (AtEvent *)eA->At(0) : nullptr;
      auto *evB = eB->GetEntriesFast() ? (AtEvent *)eB->At(0) : nullptr;
      const size_t nA = evA ? evA->GetHits().size() : 0;
      const size_t nB = evB ? evB->GetHits().size() : 0;
      nHitA += nA;
      nHitB += nB;
      if (nA != nB) {
         ++nCountMismatch;
         if (reported < maxReport) {
            printf("  event %lld: hit count %zu vs %zu\n", i, nA, nB);
            ++reported;
         }
         continue; // positions are not comparable once the counts differ
      }
      for (size_t h = 0; h < nA; ++h) {
         const auto &pA = evA->GetHits()[h]->GetPosition();
         const auto &pB = evB->GetHits()[h]->GetPosition();
         const double d = std::max({std::fabs(pA.X() - pB.X()), std::fabs(pA.Y() - pB.Y()),
                                    std::fabs(pA.Z() - pB.Z())});
         const double dq = std::fabs(evA->GetHits()[h]->GetCharge() - evB->GetHits()[h]->GetCharge());
         worst = std::max(worst, std::max(d, dq));
         if (d > tol || dq > tol) {
            ++nPosMismatch;
            if (reported < maxReport) {
               printf("  event %lld hit %zu: (%.4f,%.4f,%.4f) q=%.3f  vs  (%.4f,%.4f,%.4f) q=%.3f\n", i, h, pA.X(),
                      pA.Y(), pA.Z(), evA->GetHits()[h]->GetCharge(), pB.X(), pB.Y(), pB.Z(),
                      evB->GetHits()[h]->GetCharge());
               ++reported;
            }
         }
      }
   }

   printf("\n  events compared      : %lld\n", N);
   printf("  hits  reference/new  : %ld / %ld\n", nHitA, nHitB);
   printf("  events with differing hit count : %ld\n", nCountMismatch);
   printf("  hits differing in position/charge: %ld  (largest deviation %.3g)\n", nPosMismatch, worst);

   if (NA == NB && nCountMismatch == 0 && nPosMismatch == 0)
      printf("\n  [IDENTICAL] normal-mode output is unchanged, hit for hit.\n\n");
   else
      printf("\n  [DIFFER] the two reconstructions are NOT the same. If no deliberate change to the\n"
             "  normal path was intended, this is a regression -- check the guard in\n"
             "  AtPSA::CalculateZGeo and the ternary in AtClusterize::getCurrentPointLocation.\n\n");
}
