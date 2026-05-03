/// @file pid_confusion.C
/// @brief PID confusion matrix: MC-truth species × {UKF chosen hypothesis,
///        PRA cross-product sign}.
///
/// For each fitted track in `output_ukf_multi_riemann_1k.root`:
///   1. Pull the source AtTrack from AtTrackingEvent.GetTrackArray().
///   2. Truth species: match AtTrack hits to AtMCPoint (sim file) in (x,y)
///      within 3 mm; majority-vote MC trackID; PDG via MCTrack.
///   3. PRA cross-product sign: r²-sort hits, sign of (p_mid - p_first) ×
///      (p_last - p_first)·ẑ — calibrated convention crossZ<0 ⇒ q>0.
///   4. UKF hypothesis: AtFittedTrack.GetParticleInfo().{idPDG, charge}.
///
/// Output: ASCII confusion tables + reports per-truth-class fractions.
///
/// Run: root -b -q 'pid_confusion.C("./data/output_ukf_multi_riemann_1k.root", "./data/attpcsim.root")'

namespace {

constexpr double kMatchTol_mm = 3.0;

const std::vector<int> kTruthPDGs = {321, -321, 211, -211, 2212, 0}; // 0 = "other"
const std::vector<const char *> kTruthNames = {"K+", "K-", "pi+", "pi-", "p", "other"};

const std::vector<TString> kUKFNames = {"K+", "K-", "pi+", "pi-", "?"}; // "?" = unmatched

int truthIdx(int pdg)
{
   for (size_t i = 0; i < kTruthPDGs.size() - 1; ++i)
      if (kTruthPDGs[i] == pdg) return (int)i;
   return (int)kTruthPDGs.size() - 1;
}

int ukfIdx(const TString &name)
{
   for (size_t i = 0; i < kUKFNames.size() - 1; ++i)
      if (name == kUKFNames[i]) return (int)i;
   return (int)kUKFNames.size() - 1;
}

} // namespace

void pid_confusion(TString ukfFile = "./data/output_ukf_multi_riemann_1k.root",
                   TString simFile = "./data/attpcsim.root")
{
   auto *fU = TFile::Open(ukfFile);
   auto *fS = TFile::Open(simFile);
   if (!fU || fU->IsZombie() || !fS || fS->IsZombie()) {
      std::cerr << "Cannot open inputs\n";
      return;
   }
   auto *tU = (TTree *)fU->Get("cbmsim");
   auto *tS = (TTree *)fS->Get("cbmsim");
   if (!tU || !tS) {
      std::cerr << "Missing cbmsim trees\n";
      return;
   }

   TClonesArray *teArr = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tU->SetBranchAddress("AtTrackingEvent", &teArr);
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   const Long64_t N = std::min(tU->GetEntries(), tS->GetEntries());

   const int nTruth = (int)kTruthPDGs.size();
   const int nUKF = (int)kUKFNames.size();
   std::vector<std::vector<int>> conf(nTruth, std::vector<int>(nUKF, 0));
   std::vector<std::vector<int>> signTruth_vs_PRA(nTruth, std::vector<int>(2, 0)); // 0=q+,1=q-
   std::vector<std::vector<int>> signTruth_vs_UKF(nTruth, std::vector<int>(2, 0));

   int totalFits = 0, missingTID = 0, noTruth = 0;

   for (Long64_t e = 0; e < N; ++e) {
      tU->GetEntry(e);
      if (!teArr || teArr->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      const auto &fits = te->GetFittedTracks();
      if (fits.empty()) continue;
      auto trackArr = te->GetTrackArray();
      tS->GetEntry(e);

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
         totalFits++;
         int tid = ft->GetTrackID();
         if (tid < 0 || tid >= (int)trackArr.size()) {
            missingTID++;
            continue;
         }
         const auto &hits = trackArr[tid].GetHitArray();
         if (hits.size() < 3) continue;

         // ---- Truth PDG via (x,y) majority vote
         std::map<int, int> tidVotes;
         const double tol2 = kMatchTol_mm * kMatchTol_mm;
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
            if (bestTid >= 0) tidVotes[bestTid]++;
         }
         int truthPdg = 0;
         if (!tidVotes.empty()) {
            int bestTid = -1, bestVotes = 0;
            for (auto &kv : tidVotes)
               if (kv.second > bestVotes) { bestVotes = kv.second; bestTid = kv.first; }
            if (bestTid >= 0 && bestTid < mcTrks->GetEntries())
               truthPdg = ((AtMCTrack *)mcTrks->At(bestTid))->GetPdgCode();
         }
         if (truthPdg == 0 && tidVotes.empty()) noTruth++;
         int tIdx = truthIdx(truthPdg);

         // ---- PRA cross-product sign
         int praSign = +1;
         {
            std::vector<std::pair<double, std::pair<double, double>>> byR;
            byR.reserve(hits.size());
            for (const auto &h : hits) {
               const auto &p = h->GetPosition();
               byR.emplace_back(p.X() * p.X() + p.Y() * p.Y(),
                                std::make_pair(p.X(), p.Y()));
            }
            std::sort(byR.begin(), byR.end());
            const auto &p0 = byR.front().second;
            const auto &pm = byR[byR.size() / 2].second;
            const auto &pN = byR.back().second;
            double crossZ = (pm.first - p0.first) * (pN.second - p0.second)
                           - (pm.second - p0.second) * (pN.first - p0.first);
            praSign = (crossZ > 0) ? -1 : +1; // truth-calibrated convention
         }

         // ---- UKF chosen hypothesis
         const auto &pinfo = ft->GetParticleInfo();
         int uIdx = ukfIdx(pinfo.idPDG);
         int ukfSign = (pinfo.charge > 0) ? +1 : -1;

         conf[tIdx][uIdx]++;
         signTruth_vs_PRA[tIdx][praSign > 0 ? 0 : 1]++;
         signTruth_vs_UKF[tIdx][ukfSign > 0 ? 0 : 1]++;
      }
   }

   // ---- Print results
   std::cout << "\nTotal fits: " << totalFits << "  (missing TID: " << missingTID
             << ", no truth match: " << noTruth << ")\n";

   std::cout << "\n=== UKF hypothesis (cols) vs truth species (rows) ===\n";
   std::cout << std::setw(8) << " ";
   for (int u = 0; u < nUKF; ++u) std::cout << std::setw(8) << kUKFNames[u].Data();
   std::cout << std::setw(8) << "total\n";
   for (int t = 0; t < nTruth; ++t) {
      std::cout << std::setw(8) << kTruthNames[t];
      int tot = 0;
      for (int u = 0; u < nUKF; ++u) { std::cout << std::setw(8) << conf[t][u]; tot += conf[t][u]; }
      std::cout << std::setw(8) << tot << "\n";
   }

   std::cout << "\n=== Charge sign: truth (rows) vs PRA cross-product (cols) ===\n";
   std::cout << std::setw(8) << " " << std::setw(8) << "q+" << std::setw(8) << "q-"
             << std::setw(8) << "%right\n";
   for (int t = 0; t < nTruth; ++t) {
      int qPlus = signTruth_vs_PRA[t][0], qMinus = signTruth_vs_PRA[t][1];
      int truthQ = (kTruthPDGs[t] > 0) ? +1 : (kTruthPDGs[t] < 0 ? -1 : 0);
      double rightFrac = 0;
      if ((qPlus + qMinus) > 0)
         rightFrac = 100.0 * (truthQ > 0 ? qPlus : (truthQ < 0 ? qMinus : 0)) / (qPlus + qMinus);
      std::cout << std::setw(8) << kTruthNames[t]
                << std::setw(8) << qPlus << std::setw(8) << qMinus
                << std::setw(7) << std::fixed << std::setprecision(1) << rightFrac << "%\n";
   }

   std::cout << "\n=== Charge sign: truth (rows) vs UKF chosen (cols) ===\n";
   std::cout << std::setw(8) << " " << std::setw(8) << "q+" << std::setw(8) << "q-"
             << std::setw(8) << "%right\n";
   for (int t = 0; t < nTruth; ++t) {
      int qPlus = signTruth_vs_UKF[t][0], qMinus = signTruth_vs_UKF[t][1];
      int truthQ = (kTruthPDGs[t] > 0) ? +1 : (kTruthPDGs[t] < 0 ? -1 : 0);
      double rightFrac = 0;
      if ((qPlus + qMinus) > 0)
         rightFrac = 100.0 * (truthQ > 0 ? qPlus : (truthQ < 0 ? qMinus : 0)) / (qPlus + qMinus);
      std::cout << std::setw(8) << kTruthNames[t]
                << std::setw(8) << qPlus << std::setw(8) << qMinus
                << std::setw(7) << std::fixed << std::setprecision(1) << rightFrac << "%\n";
   }
}
