/// @file cmpReco_a1975.C
/// @brief Compare the AtPatternEvent (PRA) output of two _reco.root files event-by-event
/// — the migrated-pipeline output vs the previous-analysis production file. Reports per
/// event whether the track count + cluster count + position checksum + radius sum match.
///   root -l -b -q 'cmpReco_a1975.C("run_0106", "/mnt/f/a1975/reco_test/", "/mnt/f/a1975/reco/", 3000)'

void sig(TTree *t, std::vector<std::array<double, 4>> &out, Long64_t nEv)
{
   auto *arr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &arr);
   Long64_t n = std::min(nEv, t->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)arr->At(0);
      double nt = 0, nc = 0, cks = 0, rs = 0;
      if (pe)
         for (auto &tr : pe->GetTrackCand()) {
            nt += 1;
            auto *cl = tr.GetHitClusterArray();
            nc += cl->size();
            for (auto &c : *cl) { auto p = c.GetPosition(); cks += p.X() + p.Y() + p.Z(); }
            rs += tr.GetGeoRadius();
         }
      out.push_back({nt, nc, cks, rs});
   }
}

void cmpReco_a1975(TString runName = "run_0106", TString newDir = "/mnt/f/a1975/reco_test/",
                   TString prevDir = "/mnt/f/a1975/reco/", Long64_t nEv = 3000)
{
   gSystem->Load("libAtReconstruction.so");
   TFile *fa = TFile::Open(newDir + runName + "_reco.root");
   TFile *fb = TFile::Open(prevDir + runName + "_reco.root");
   if (!fa || !fb) { std::cout << runName << ": missing file(s)\n"; return; }
   std::vector<std::array<double, 4>> A, B;
   sig((TTree *)fa->Get("cbmsim"), A, nEv);
   sig((TTree *)fb->Get("cbmsim"), B, nEv);
   Long64_t n = std::min(A.size(), B.size());
   int diffEv = 0, trA = 0, trB = 0;
   double maxRelRad = 0;
   for (Long64_t i = 0; i < n; ++i) {
      trA += (int)A[i][0]; trB += (int)B[i][0];
      bool d = A[i][0] != B[i][0] || A[i][1] != B[i][1] || std::abs(A[i][2] - B[i][2]) > 1e-3;
      double rd = (B[i][3] != 0) ? std::abs(A[i][3] - B[i][3]) / std::abs(B[i][3]) : 0;
      maxRelRad = std::max(maxRelRad, rd);
      if (d) { diffEv++; if (diffEv <= 6) printf("  DIFF ev %lld: new(t=%g,c=%g) prev(t=%g,c=%g)\n", i, A[i][0], A[i][1], B[i][0], B[i][1]); }
   }
   printf("%s: %lld events | tracks new=%d prev=%d | track/cluster/pos mismatches=%d | max rel radius diff=%.2e %s\n",
          runName.Data(), n, trA, trB, diffEv, maxRelRad,
          (diffEv == 0 ? "\033[1;32m[MATCH]\033[0m" : "\033[1;31m[DIFFERS]\033[0m"));
}
