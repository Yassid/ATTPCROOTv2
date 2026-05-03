/// @file track_loss.C
/// @brief PRA → UKF track yield diagnostic.
///
/// Run: root -b -q track_loss.C

void track_loss(TString digiFile = "./data/output_digi.root",
                TString ukfFile = "./data/output_ukf_multi.root")
{
   const int NB = 13;
   std::vector<int> binPRA(NB, 0);
   int totPRA = 0;
   int evPRA = 0;

   // Per-event PRA hit counts so we can later look up by GetTrackID() when
   // counting fits. (Matches across files since UKF runs on output_digi.)
   std::vector<std::vector<int>> hitCountByEvent;

   // Open the sim file and pre-attach AtMCPoint/AtMCTrack — without this,
   // Cling crashes when iterating AtPatternEvent.GetTrackCand().GetHitArray()
   // (presumably AtData dictionary load order).
   auto *fSim = TFile::Open("./data/attpcsim.root");
   auto *tSim = (TTree *)fSim->Get("cbmsim");
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tSim->SetBranchAddress("MCTrack", &mcTrks);

   // ---- Pass 1: digi/PRA only
   auto *fDigi = TFile::Open(digiFile);
   auto *tDigi = (TTree *)fDigi->Get("cbmsim");
   TClonesArray *peArr = nullptr;
   tDigi->SetBranchAddress("AtPatternEvent", &peArr);
   const Long64_t nDigi = tDigi->GetEntries();
   hitCountByEvent.resize(nDigi);
   for (Long64_t e = 0; e < nDigi; ++e) {
      tDigi->GetEntry(e);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      auto &tracks = pe->GetTrackCand();
      if (!tracks.empty()) evPRA++;
      hitCountByEvent[e].reserve(tracks.size());
      for (auto &tr : tracks) {
         int n = (int)tr.GetHitArray().size();
         hitCountByEvent[e].push_back(n);
         int b = std::min(n / 5, NB - 1);
         binPRA[b]++;
         totPRA++;
      }
   }

   // ---- Pass 2: UKF only
   auto *fUkf = TFile::Open(ukfFile);
   auto *tUkf = (TTree *)fUkf->Get("cbmsim");
   TClonesArray *teArr = nullptr;
   tUkf->SetBranchAddress("AtTrackingEvent", &teArr);
   const Long64_t nUkf = tUkf->GetEntries();
   const Long64_t N = std::min(nDigi, nUkf);

   std::vector<int> binFit(NB, 0);
   int totFit = 0, evFit = 0;
   for (Long64_t e = 0; e < N; ++e) {
      tUkf->GetEntry(e);
      auto *te = (AtTrackingEvent *)teArr->At(0);
      const auto &fits = te->GetFittedTracks();
      if (!fits.empty()) evFit++;
      totFit += (int)fits.size();
      for (auto &ft : fits) {
         int tid = ft->GetTrackID();
         if (tid >= 0 && tid < (int)hitCountByEvent[e].size()) {
            int n = hitCountByEvent[e][tid];
            int b = std::min(n / 5, NB - 1);
            binFit[b]++;
         }
      }
   }

   std::cout << "\nevents processed   = " << N << "\n";
   std::cout << "events with PRA    = " << evPRA << " (" << 100.0 * evPRA / N << "%)\n";
   std::cout << "events with UKF    = " << evFit << " (" << 100.0 * evFit / N << "%)\n";
   std::cout << "PRA tracks total   = " << totPRA << "\n";
   std::cout << "UKF fits total     = " << totFit
             << "  (" << 100.0 * totFit / std::max(totPRA, 1) << "% of PRA)\n";

   std::cout << "\nAcceptance by PRA hit count:\n";
   std::cout << "   bin     PRA     fit  fit/pra\n";
   for (int b = 0; b < NB; ++b) {
      if (binPRA[b] == 0 && binFit[b] == 0) continue;
      double f = (binPRA[b] > 0) ? (double)binFit[b] / binPRA[b] : 0.0;
      int lo = b * 5;
      int hi = (b == NB - 1) ? 9999 : (lo + 4);
      char hiBuf[8];
      if (hi == 9999) snprintf(hiBuf, 8, "INF");
      else snprintf(hiBuf, 8, "%d", hi);
      std::cout << "   " << std::setw(2) << lo << "-" << std::setw(3) << hiBuf
                << "  " << std::setw(6) << binPRA[b]
                << "  " << std::setw(6) << binFit[b]
                << "   " << std::fixed << std::setprecision(2) << f << "\n";
   }
}
