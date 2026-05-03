/// @file theta_bias.C
/// @brief Compare smoothed θ at first cluster (Kinematics.theta) to MC truth
///        θ for every truth-matched fit. Reports the (Δθ, Δφ) distribution
///        and the per-bin Δz_vertex correlation, to test whether the +8 mm
///        Δz bias comes from a θ shift toward 90°.
///
/// Run: root -b -q theta_bias.C
///      root -b -q 'theta_bias.C("./data/output_digi.root", "./data/attpcsim.root")'

void theta_bias(TString digiFile = "./data/output_digi.root",
                TString simFile = "./data/attpcsim.root")
{
   auto *fD = TFile::Open(digiFile);
   auto *fS = TFile::Open(simFile);
   if (!fD || fD->IsZombie() || !fS || fS->IsZombie()) {
      std::cerr << "Cannot open inputs\n"; return;
   }
   auto *tD = (TTree *)fD->Get("cbmsim");
   auto *tS = (TTree *)fS->Get("cbmsim");

   TClonesArray *teA = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tD->SetBranchAddress("AtTrackingEvent", &teA);
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   const Long64_t N = std::min(tD->GetEntries(), tS->GetEntries());

   std::vector<double> dTheta, dPhi, dZ, truthTheta;
   int matched = 0, unmatched = 0;

   for (Long64_t e = 0; e < N; ++e) {
      tD->GetEntry(e);
      tS->GetEntry(e);
      if (!teA || teA->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teA->At(0);
      const auto &fits = te->GetFittedTracks();
      if (fits.empty()) continue;
      auto trArr = te->GetTrackArray();

      // Cache MC points (mm) for (x,y)-match
      const int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC);
      std::vector<int> mcTid(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10.0;
         mcY[k] = mp->GetY() * 10.0;
         mcTid[k] = mp->GetTrackID();
      }

      for (const auto &ft : fits) {
         int tid = ft->GetTrackID();
         if (tid < 0 || tid >= (int)trArr.size()) continue;
         const auto &hits = trArr[tid].GetHitArray();
         if (hits.empty()) continue;

         // (x,y) majority-vote → MC primary tid
         std::map<int, int> votes;
         const double tol2 = 9.0; // 3 mm
         for (const auto &h : hits) {
            const auto &p = h->GetPosition();
            double bestD2 = tol2;
            int bestTid = -1;
            for (int k = 0; k < nMC; ++k) {
               double dx = p.X() - mcX[k];
               double dy = p.Y() - mcY[k];
               double d2 = dx * dx + dy * dy;
               if (d2 < bestD2) { bestD2 = d2; bestTid = mcTid[k]; }
            }
            if (bestTid >= 0) votes[bestTid]++;
         }
         int truthTid = -1, bestN = 0;
         for (auto &kv : votes) if (kv.second > bestN) { bestN = kv.second; truthTid = kv.first; }
         if (truthTid < 0 || truthTid >= mcTrks->GetEntries()) { unmatched++; continue; }
         matched++;

         auto *mt = (AtMCTrack *)mcTrks->At(truthTid);
         double pmag = std::sqrt(mt->GetPx() * mt->GetPx() + mt->GetPy() * mt->GetPy()
                                 + mt->GetPz() * mt->GetPz());
         if (pmag <= 0) continue;
         double thT = std::acos(mt->GetPz() / pmag) * 180.0 / TMath::Pi();
         double phT = std::atan2(mt->GetPy(), mt->GetPx()) * 180.0 / TMath::Pi();

         const auto &kin = ft->GetKinematics();
         double thF = kin.theta * 180.0 / TMath::Pi();
         double phF = kin.phi * 180.0 / TMath::Pi();
         double dphi = phF - phT;
         while (dphi > 180) dphi -= 360;
         while (dphi < -180) dphi += 360;

         const auto &props = ft->GetTrackPropertiesStruct();
         double zVtxFit = props.initialPositionXtr.Z();
         double zVtxMC = mt->GetStartZ() * 10.0;
         double dz = zVtxFit - zVtxMC;

         dTheta.push_back(thF - thT);
         dPhi.push_back(dphi);
         dZ.push_back(dz);
         truthTheta.push_back(thT);
      }
   }

   auto med = [](std::vector<double> v) {
      if (v.empty()) return 0.0;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };
   auto sigmaIQR = [](std::vector<double> v) {
      if (v.size() < 4) return 0.0;
      std::sort(v.begin(), v.end());
      return (v[(3 * v.size()) / 4] - v[v.size() / 4]) / 1.349;
   };

   std::cout << "\nMatched fits: " << matched << ", unmatched: " << unmatched << "\n";
   std::cout << "Δθ (fit - truth) [deg]:  median = " << med(dTheta)
             << ",  σ_IQR = " << sigmaIQR(dTheta) << "\n";
   std::cout << "Δφ (fit - truth) [deg]:  median = " << med(dPhi)
             << ",  σ_IQR = " << sigmaIQR(dPhi) << "\n";
   std::cout << "Δz (fit - truth) [mm]:   median = " << med(dZ)
             << ",  σ_IQR = " << sigmaIQR(dZ) << "\n";

   // Bin Δθ and Δz by truth θ to see directionality.
   std::cout << "\n  truth θ band      n     <Δθ>    σΔθ      <Δz>    σΔz\n";
   const std::vector<std::pair<double, double>> bands = {
      {30, 60}, {60, 75}, {75, 90}, {90, 105}, {105, 120}, {120, 150}};
   for (auto &b : bands) {
      std::vector<double> th, zz;
      for (size_t i = 0; i < truthTheta.size(); ++i) {
         if (truthTheta[i] >= b.first && truthTheta[i] < b.second) {
            th.push_back(dTheta[i]);
            zz.push_back(dZ[i]);
         }
      }
      std::cout << "  [" << std::setw(3) << (int)b.first << ","
                << std::setw(3) << (int)b.second << ")  "
                << std::setw(5) << th.size() << "  "
                << std::fixed << std::setprecision(2)
                << std::setw(7) << med(th) << "  "
                << std::setw(6) << sigmaIQR(th) << "   "
                << std::setw(7) << med(zz) << "  "
                << std::setw(6) << sigmaIQR(zz) << "\n";
   }
}
