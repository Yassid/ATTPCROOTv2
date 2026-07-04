/// @file charge_test.C
/// @brief Offline A/B of charge-sign estimators on the branch-8 PRA tracks, vs truth.
///   OLD  = 3-point cross product (nearest/mid/farthest by r) — the current AtPRA rule.
///   NEW  = angular sweep around the FITTED circle centre (GetGeoCenter): robust because
///          |p - C| ~ R (~312 mm) gives a large lever arm on shallow arcs.
/// No rebuild needed — both are computed from the stored hits + GeoCenter.
/// Run: root -b -q charge_test.C
void charge_test(TString digiFile = "./data/output_digi_both8.root", TString simFile = "./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   const double kTol = 3.0;
   TFile fD(digiFile); TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);  TTree *tS = (TTree *)fS.Get("cbmsim");
   TClonesArray *pat = new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent", &pat);
   TClonesArray *mcPts = new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint", &mcPts);
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack"); tS->SetBranchAddress("MCTrack", &mcTrks);

   // accuracy counters for each estimator + sign convention
   long tot = 0;
   long oldRight = 0, newRight = 0;
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e); tS->GetEntry(e);
      if (!pat->GetEntries()) continue;
      auto *pe = (AtPatternEvent *)pat->At(0);
      int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC); std::vector<int> mcPdg(nMC);
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k); mcX[k] = mp->GetX() * 10; mcY[k] = mp->GetY() * 10;
         int t = mp->GetTrackID(); auto *mt = (t >= 0 && t < mcTrks->GetEntries()) ? (AtMCTrack *)mcTrks->At(t) : nullptr; mcPdg[k] = mt ? mt->GetPdgCode() : 0; }

      for (auto &tr : pe->GetTrackCand()) {
         const auto &hits = tr.GetHitArray();
         if (hits.size() < 3) continue;
         // order hits by distance from origin (vertex first)
         std::vector<std::pair<double, std::pair<double, double>>> byR;
         for (const auto &h : hits) { const auto &p = h->GetPosition(); byR.emplace_back(p.X() * p.X() + p.Y() * p.Y(), std::make_pair(p.X(), p.Y())); }
         std::sort(byR.begin(), byR.end());
         const auto &p0 = byR.front().second; const auto &pm = byR[byR.size() / 2].second; const auto &pN = byR.back().second;

         // OLD: 3-point cross (current convention: crossZ<0 => q>0)
         double crossOld = (pm.first - p0.first) * (pN.second - p0.second) - (pm.second - p0.second) * (pN.first - p0.first);
         int qOld = (crossOld > 0) ? -1 : +1;

         // NEW: sweep around fitted centre
         auto c = tr.GetGeoCenter();
         double v0x = p0.first - c.first, v0y = p0.second - c.second;
         double vNx = pN.first - c.first, vNy = pN.second - c.second;
         double crossNew = v0x * vNy - v0y * vNx;
         int qNew = (crossNew > 0) ? -1 : +1; // provisional same convention; calibrated below

         // truth charge via (x,y) hit majority vote
         std::map<int, int> votes;
         for (const auto &h : hits) { const auto &p = h->GetPosition(); double best = kTol * kTol; int bp = 0;
            for (int k = 0; k < nMC; ++k) { double d2 = (p.X() - mcX[k]) * (p.X() - mcX[k]) + (p.Y() - mcY[k]) * (p.Y() - mcY[k]); if (d2 < best) { best = d2; bp = mcPdg[k]; } }
            if (bp) votes[bp]++; }
         int tp = 0, bv = 0; for (auto &kv : votes) if (kv.second > bv) { bv = kv.second; tp = kv.first; }
         if (!tp) continue;
         int qTruth = (tp > 0) ? +1 : -1; // pi+ = q>0
         tot++;
         if (qOld == qTruth) oldRight++;
         if (qNew == qTruth) newRight++;
      }
   }
   double accOld = 100.0 * oldRight / tot;
   double accNew = 100.0 * newRight / tot;
   // if NEW came out anti-correlated, the convention is flipped — report the flipped accuracy
   printf("tracks with truth: %ld\n", tot);
   printf("OLD (3-point)          : %.1f%%\n", accOld);
   printf("NEW (sweep-about-centre): %.1f%%%s\n", std::max(accNew, 100 - accNew),
          (accNew < 50) ? "   [sign convention flipped vs OLD]" : "");
}
