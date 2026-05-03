/// @file digi_z_offset.C
/// @brief Compare MC truth point z to reconstructed digi hit z, matched by
///        (x,y) proximity. Reports the population mean/median Δz between
///        digi and MC, to test whether the +8 mm Δz vertex bias is rooted
///        in the digi/PSA z reconstruction (not the back-extrap).
///
/// Run: root -b -q digi_z_offset.C

void digi_z_offset(TString digiFile = "./data/output_digi.root",
                   TString simFile = "./data/attpcsim.root")
{
   auto *fD = TFile::Open(digiFile);
   auto *fS = TFile::Open(simFile);
   if (!fD || fD->IsZombie() || !fS || fS->IsZombie()) {
      std::cerr << "Cannot open inputs\n"; return;
   }
   auto *tD = (TTree *)fD->Get("cbmsim");
   auto *tS = (TTree *)fS->Get("cbmsim");

   TClonesArray *peA = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   tD->SetBranchAddress("AtPatternEvent", &peA);
   tS->SetBranchAddress("AtTpcPoint", &mcPts);

   const Long64_t N = std::min(tD->GetEntries(), tS->GetEntries());

   std::vector<double> dZ;
   std::vector<double> zMC, zDigi;
   const double tol2 = 4.0; // 2 mm xy match tolerance

   for (Long64_t e = 0; e < N; ++e) {
      tD->GetEntry(e);
      tS->GetEntry(e);
      if (!peA || peA->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)peA->At(0);
      const auto &tracks = pe->GetTrackCand();
      if (tracks.empty()) continue;

      // Cache MC points (mm)
      const int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC), mcZ(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10.0;
         mcY[k] = mp->GetY() * 10.0;
         mcZ[k] = mp->GetZ() * 10.0;
      }

      // For each digi hit, find the closest MC point in xy
      for (const auto &tr : tracks) {
         for (auto &hp : tr.GetHitArray()) {
            const auto &p = hp->GetPosition();
            double bestD2 = tol2;
            int bestK = -1;
            for (int k = 0; k < nMC; ++k) {
               double dx = p.X() - mcX[k];
               double dy = p.Y() - mcY[k];
               double d2 = dx * dx + dy * dy;
               if (d2 < bestD2) { bestD2 = d2; bestK = k; }
            }
            if (bestK >= 0) {
               dZ.push_back(p.Z() - mcZ[bestK]);
               zDigi.push_back(p.Z());
               zMC.push_back(mcZ[bestK]);
            }
         }
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

   std::cout << "\nMatched hit↔MC pairs: " << dZ.size() << "\n";
   std::cout << "Δz (digi - MC) [mm]:  median = " << med(dZ)
             << ",  σ_IQR = " << sigmaIQR(dZ) << "\n";

   // Bin Δz by MC z to see whether the offset is constant (additive) or
   // depends on z (multiplicative drift-velocity error).
   std::cout << "\n   MC z band         n     <Δz>     σΔz\n";
   const std::vector<std::pair<double, double>> bands = {
      {0, 30}, {30, 60}, {60, 90}, {90, 120}, {120, 150}, {150, 180}};
   for (auto &b : bands) {
      std::vector<double> v;
      for (size_t i = 0; i < dZ.size(); ++i)
         if (zMC[i] >= b.first && zMC[i] < b.second) v.push_back(dZ[i]);
      std::cout << "  [" << std::setw(3) << (int)b.first << ","
                << std::setw(3) << (int)b.second << ")  "
                << std::setw(6) << v.size() << "  "
                << std::fixed << std::setprecision(2)
                << std::setw(7) << med(v) << "   "
                << std::setw(6) << sigmaIQR(v) << "\n";
   }
}
