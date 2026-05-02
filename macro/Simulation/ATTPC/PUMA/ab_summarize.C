/// @file ab_summarize.C
/// @brief Summarize one PRA/UKF run for the Riemann-vs-TC A/B comparison.
///
/// Reads ./data/output_digi<tag>.root and ./data/output_ukf_multi<tag>.root and
/// prints PRA candidates, UKF fits, and hypothesis distribution.
///
/// Run: root -b -q 'ab_summarize.C("_TC")'

void ab_summarize(TString tag = "")
{
   TString digiFile = TString("./data/output_digi") + tag + ".root";
   TString ukfFile = TString("./data/output_ukf_multi") + tag + ".root";

   // --- PRA stats ---
   TFile *fPra = TFile::Open(digiFile);
   if (!fPra || fPra->IsZombie()) {
      std::cerr << "Cannot open " << digiFile << "\n";
      return;
   }
   TTree *tPra = dynamic_cast<TTree *>(fPra->Get("cbmsim"));
   TClonesArray *patArr = nullptr;
   tPra->SetBranchAddress("AtPatternEvent", &patArr);

   Long64_t nev = tPra->GetEntries();
   int eventsWithCand = 0;
   int totalCand = 0;
   for (Long64_t i = 0; i < nev; ++i) {
      tPra->GetEntry(i);
      auto *patEvt = dynamic_cast<AtPatternEvent *>(patArr->At(0));
      if (!patEvt)
         continue;
      const auto &tracks = patEvt->GetTrackCand();
      int n = tracks.size();
      if (n > 0)
         ++eventsWithCand;
      totalCand += n;
   }
   fPra->Close();

   // --- UKF stats ---
   TFile *fUkf = TFile::Open(ukfFile);
   if (!fUkf || fUkf->IsZombie()) {
      std::cerr << "Cannot open " << ukfFile << "\n";
      return;
   }
   TTree *tUkf = dynamic_cast<TTree *>(fUkf->Get("cbmsim"));
   TClonesArray *trkArr = nullptr;
   tUkf->SetBranchAddress("AtTrackingEvent", &trkArr);

   Long64_t nevU = tUkf->GetEntries();
   int eventsWithFit = 0;
   int totalFits = 0;
   std::map<std::string, int> hypoCount;
   for (Long64_t i = 0; i < nevU; ++i) {
      tUkf->GetEntry(i);
      auto *trkEvt = dynamic_cast<AtTrackingEvent *>(trkArr->At(0));
      if (!trkEvt)
         continue;
      const auto &fits = trkEvt->GetFittedTracks();
      int n = fits.size();
      if (n > 0)
         ++eventsWithFit;
      totalFits += n;
      for (const auto &ft : fits) {
         std::string name = std::string(ft.GetParticleInfo().idPDG.Data());
         if (name.empty())
            name = "<unset>";
         ++hypoCount[name];
      }
   }
   fUkf->Close();

   std::cout << "\n=== A/B summary tag='" << tag << "' ===\n";
   std::cout << "  events (digi/ukf)        : " << nev << " / " << nevU << "\n";
   std::cout << "  events with PRA cand.    : " << eventsWithCand << " (" << totalCand << " total cands)\n";
   std::cout << "  events with UKF fit      : " << eventsWithFit << " (" << totalFits << " total fits)\n";
   std::cout << "  hypothesis distribution  :\n";
   for (auto &kv : hypoCount) {
      std::cout << "      " << std::setw(8) << std::left << kv.first << " : " << kv.second << "\n";
   }
   std::cout << "  truth (per event)        : 2 K+ + 1 pi-\n";
}
