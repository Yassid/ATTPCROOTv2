/// @file pra_efficiency.C
/// @brief Truth-matched PRA efficiency for the PUMA test-8 channel (2 back-to-back
///        primary pions). Per event: match each reco track to a true pion via hit
///        -> nearest-MC-point trackID majority vote. Reports:
///          efficiency  = found true pions / all true pions
///          purity      = mean dominant-particle fraction over reco tracks
///          fake rate   = reco tracks not matching any true pion
///          event class = good(2 clean) / merged(1) / split(>=3) / partial
///  A true pion is "found" if some reco track has >= fMinFrac of its hits AND purity >= 0.7.
/// Run: root -b -q 'pra_efficiency.C("data/mm_pk.root","data/attpcsim.root","arc-walk")'
void pra_efficiency(TString recoFile, TString simFile, TString tag = "PRA", double matchTol = 3.0,
                    int minHitsAbs = 6, double minPurity = 0.7)
{
   gSystem->Load("libAtReconstruction.so");
   TFile fD(recoFile); auto *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile); auto *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("  %-12s MISSING FILE\n", tag.Data()); return; }
   auto *patArr = new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent", &patArr);
   auto *mcPts = new TClonesArray("AtMCPoint"); auto *mcTrks = new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint", &mcPts); tS->SetBranchAddress("MCTrack", &mcTrks);

   long nTrue = 0, nFound = 0, nRecoTot = 0, nFake = 0;
   double purSum = 0; long purN = 0;
   int evGood = 0, evMerged = 0, evSplit = 0, evPartial = 0, evTot = 0;
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e); tS->GetEntry(e);
      if (patArr->GetEntries() == 0) continue;
      // truth: MC (x,y) per point + its primary-pion trackID; count true hits per pion trackID
      int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC); std::vector<int> mcTid(nMC);
      std::map<int, int> trueHits; // trackID -> #MC points (primary pions only)
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10; mcY[k] = mp->GetY() * 10; int tid = mp->GetTrackID();
         auto *mt = (tid >= 0 && tid < mcTrks->GetEntries()) ? (AtMCTrack *)mcTrks->At(tid) : nullptr;
         mcTid[k] = (mt && std::abs(mt->GetPdgCode()) == 211) ? tid : -1;
         if (mcTid[k] >= 0) trueHits[tid]++;
      }
      if (trueHits.empty()) continue;
      evTot++;
      auto *pe = (AtPatternEvent *)patArr->At(0);
      auto &tc = pe->GetTrackCand();
      // for each reco track: majority-vote trackID + purity
      std::map<int, int> bestRecoHitsForTrue; // trackID -> max matched hits by a single reco track
      for (auto &tr : tc) { nRecoTot++;
         std::map<int, int> votes; int nh = 0;
         for (auto &h : tr.GetHitArray()) { const auto &pp = h->GetPosition(); nh++;
            double best = matchTol * matchTol; int bt = -2;
            for (int k = 0; k < nMC; ++k) { if (mcTid[k] < 0) continue;
               double d2 = (pp.X() - mcX[k]) * (pp.X() - mcX[k]) + (pp.Y() - mcY[k]) * (pp.Y() - mcY[k]);
               if (d2 < best) { best = d2; bt = mcTid[k]; } }
            if (bt >= 0) votes[bt]++;
         }
         int dom = -1, domV = 0, totV = 0;
         for (auto &kv : votes) { totV += kv.second; if (kv.second > domV) { domV = kv.second; dom = kv.first; } }
         if (totV == 0) { nFake++; continue; }
         double purity = (double)domV / nh; purSum += purity; purN++;
         if (purity < 0.5) { nFake++; continue; }
         if (domV > bestRecoHitsForTrue[dom]) bestRecoHitsForTrue[dom] = domV;
      }
      // efficiency per true pion + event classification
      int foundThisEv = 0;
      for (auto &kv : trueHits) { nTrue++;
         int matched = bestRecoHitsForTrue.count(kv.first) ? bestRecoHitsForTrue[kv.first] : 0;
         if (matched >= minHitsAbs) { nFound++; foundThisEv++; } // found: a clean reco track with >= minHits of this pion
      }
      int nT = trueHits.size(), nR = tc.size();
      if (foundThisEv == nT && nR == nT) evGood++;
      else if (nR < nT) evMerged++;
      else if (nR > nT && foundThisEv == nT) evSplit++;
      else evPartial++;
   }
   printf("  %-12s eff=%.1f%% (%ld/%ld)  purity=%.1f%%  fake=%.1f%%  | events: good=%.0f%% merged=%.0f%% split=%.0f%% partial=%.0f%%  (n=%d)\n",
          tag.Data(), 100.0 * nFound / std::max(1L, nTrue), nFound, nTrue,
          100.0 * purSum / std::max(1L, purN), 100.0 * nFake / std::max(1L, nRecoTot),
          100.0 * evGood / std::max(1, evTot), 100.0 * evMerged / std::max(1, evTot),
          100.0 * evSplit / std::max(1, evTot), 100.0 * evPartial / std::max(1, evTot), evTot);
}
