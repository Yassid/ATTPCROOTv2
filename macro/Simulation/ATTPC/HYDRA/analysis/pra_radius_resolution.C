/// @file pra_radius_resolution.C
/// @brief Measure σ(R)/R of the PRA-seed circle radius vs MC truth.
///
/// Usage:
///   root -b -q 'analysis/pra_radius_resolution.C("data/output_digi.root","baseline")'
///   root -b -q 'analysis/pra_radius_resolution.C("data/output_digi_daw2.root","drift-aware")'

void pra_radius_resolution(const char *digiFile = "data/output_digi.root",
                           const char *label = "baseline",
                           const char *simFile = "data/attpcsim.root",
                           double Bz_T = 0.5)
{
   TFile fSim(simFile);
   TFile fDigi(digiFile);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   auto *trks = new TClonesArray("AtMCTrack");
   auto *patEvtArr = new TClonesArray("AtPatternEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tDigi->SetBranchAddress("AtPatternEvent", &patEvtArr);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<TH1F *> hDR(NB);
   for (int b = 0; b < NB; ++b)
      hDR[b] = new TH1F(Form("hDR_%d", b), "", 100, -0.5, 0.5);

   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;

      double px = mc->GetPx(), py = mc->GetPy(), pz = mc->GetPz();
      double p = std::sqrt(px * px + py * py + pz * pz); // GeV
      double pT_MeV = std::sqrt(px * px + py * py) * 1000.0;
      double R_true = pT_MeV / (0.3 * Bz_T); // mm
      double thMC = std::acos(pz / p) * 180. / M_PI;

      int b = -1;
      for (int k = 0; k < NB; ++k)
         if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;

      if (patEvtArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patEvtArr->At(0);
      auto &tracks = pat->GetTrackCand();
      if (tracks.empty()) continue;
      // Take the largest track candidate
      const AtTrack *best = nullptr;
      size_t bestN = 0;
      for (const auto &t : tracks) {
         size_t nHits = t.GetHitArray().size();
         if (nHits > bestN) { bestN = nHits; best = &t; }
      }
      if (!best) continue;
      double R_fit = best->GetGeoRadius();
      if (R_fit <= 0 || !std::isfinite(R_fit)) continue;
      hDR[b]->Fill((R_fit - R_true) / R_true);
   }

   std::cout << "\n=== pra_radius_resolution [" << label << "] ===\n";
   std::cout << "theta(deg)   N       <ΔR/R>     σ(R)/R\n";
   std::cout << std::string(50, '-') << "\n";
   double allMean = 0, allSig = 0, allN = 0;
   for (int b = 0; b < NB; ++b) {
      if (hDR[b]->GetEntries() < 15) {
         printf("%4.0f-%-4.0f   %4.0f   ----        ----\n",
                edges[b], edges[b + 1], hDR[b]->GetEntries());
         continue;
      }
      double mean = hDR[b]->GetMean(), rms = hDR[b]->GetRMS();
      hDR[b]->Fit("gaus", "Q0", "", mean - 2.5 * rms, mean + 2.5 * rms);
      auto *f = hDR[b]->GetFunction("gaus");
      double mu = f ? f->GetParameter(1) : mean;
      double sig = f ? f->GetParameter(2) : rms;
      printf("%4.0f-%-4.0f   %4.0f   %+6.3f%%    %5.2f%%\n",
             edges[b], edges[b + 1], hDR[b]->GetEntries(), mu * 100., sig * 100.);
      allMean += mu * hDR[b]->GetEntries();
      allSig += sig * hDR[b]->GetEntries();
      allN += hDR[b]->GetEntries();
   }
   if (allN > 0)
      printf("ALL          %4.0f   %+6.3f%%    %5.2f%%\n", allN, allMean / allN * 100., allSig / allN * 100.);
}
