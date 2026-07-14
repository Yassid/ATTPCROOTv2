/// @file track_fit_efficiency.C
/// @brief Tracking + fitting efficiency for a PUMA reco file. Reports, per event:
///   tracking eff  = true primaries found by pattern recognition / all true primaries
///   PRA tracks/evt
///   fitting eff   = successfully fitted tracks / PRA tracks   (UKF and GENFIT)
///   overall eff   = fitted-and-truth-matched / all true primaries   (end-to-end)
/// Truth match: hit -> nearest MC point trackID majority (purity >= 0.5, >= minHits).
/// Run: root -b -q 'track_fit_efficiency.C("data/hs_reco375.root","data/hs_sim375.root","DLC 375")'
void track_fit_efficiency(TString recoFile, TString simFile, TString tag = "reco",
                          double matchTol = 3.0, int minHits = 6, int targetPDG = 211)
{
   gSystem->Load("libAtReconstruction.so");
   TFile fD(recoFile); auto *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile); auto *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("  %-14s MISSING\n", tag.Data()); return; }
   auto *patArr = new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent", &patArr);
   auto *ukfArr = new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   auto *gfArr = new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   auto *mcPts = new TClonesArray("AtMCPoint"); auto *mcTrks = new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint", &mcPts); tS->SetBranchAddress("MCTrack", &mcTrks);

   long nTrue = 0, nPRA = 0, nPRAfound = 0, nUKF = 0, nGF = 0; int evTot = 0;
   auto nFits = [](TClonesArray *a) { if (a->GetEntries() == 0) return 0; int n = 0;
      for (const auto &ft : ((AtTrackingEvent *)a->At(0))->GetFittedTracks()) if (ft->GetKinematics(0).kineticEnergy > 0) n++;
      return n; };
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < nE; ++e) { tD->GetEntry(e); tS->GetEntry(e);
      if (patArr->GetEntries() == 0) continue;
      int nMC = mcPts->GetEntries(); std::vector<double> mcX(nMC), mcY(nMC); std::vector<int> mcTid(nMC);
      std::map<int, int> trueHits;
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k); mcX[k] = mp->GetX() * 10; mcY[k] = mp->GetY() * 10;
         int tid = mp->GetTrackID(); auto *mt = (tid >= 0 && tid < mcTrks->GetEntries()) ? (AtMCTrack *)mcTrks->At(tid) : nullptr;
         mcTid[k] = (mt && std::abs(mt->GetPdgCode()) == targetPDG) ? tid : -1; if (mcTid[k] >= 0) trueHits[tid]++; }
      if (trueHits.empty()) continue; evTot++;
      auto *pe = (AtPatternEvent *)patArr->At(0); auto &tc = pe->GetTrackCand();
      std::map<int, int> bestHits;
      for (auto &tr : tc) { nPRA++; std::map<int, int> votes; int nh = 0;
         for (auto &h : tr.GetHitArray()) { const auto &pp = h->GetPosition(); nh++; double best = matchTol * matchTol; int bt = -1;
            for (int k = 0; k < nMC; ++k) { if (mcTid[k] < 0) continue; double d2 = (pp.X()-mcX[k])*(pp.X()-mcX[k])+(pp.Y()-mcY[k])*(pp.Y()-mcY[k]);
               if (d2 < best) { best = d2; bt = mcTid[k]; } } if (bt >= 0) votes[bt]++; }
         int dom = -1, dv = 0, tv = 0; for (auto &kv : votes) { tv += kv.second; if (kv.second > dv) { dv = kv.second; dom = kv.first; } }
         if (tv > 0 && (double)dv / nh >= 0.5 && dv > bestHits[dom]) bestHits[dom] = dv; }
      for (auto &kv : trueHits) { nTrue++; if (bestHits.count(kv.first) && bestHits[kv.first] >= minHits) nPRAfound++; }
      nUKF += nFits(ukfArr); nGF += nFits(gfArr);
   }
   double trk = 100.0 * nPRAfound / std::max(1L, nTrue);
   printf("  %-14s | events %d | true %ld | PRA %.2f/evt | TRACKING eff %.1f%% | FIT eff UKF %.1f%% GENFIT %.1f%% | end-to-end UKF %.1f%% GENFIT %.1f%%\n",
          tag.Data(), evTot, nTrue, (double)nPRA / evTot, trk,
          100.0 * nUKF / std::max(1L, nPRA), 100.0 * nGF / std::max(1L, nPRA),
          100.0 * std::min((long)nUKF, nPRAfound) / std::max(1L, nTrue), 100.0 * std::min((long)nGF, nPRAfound) / std::max(1L, nTrue));
}
