/// @file hit_residual.C
/// @brief Per-hit perpendicular residual to the PRA-fitted circle, binned
/// by drift z. Diagnoses whether the σ(R)/R floor is set by Gaussian per-hit
/// σ_xy (which should be ~1 mm, growing with √z due to diffusion) or by
/// correlated noise (residuals dominated by track-level shifts).
///
/// Usage:
///   root -b -q 'analysis/hit_residual.C("data/output_digi.root")'

void hit_residual(const char *digiFile = "data/output_digi.root")
{
   TFile fDigi(digiFile);
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   auto *patEvtArr = new TClonesArray("AtPatternEvent");
   tDigi->SetBranchAddress("AtPatternEvent", &patEvtArr);

   // residuals binned by z range
   const int NZ = 5;
   double zEdges[NZ + 1] = {0., 200., 400., 600., 800., 1000.};
   std::vector<TH1F *> hRes(NZ);
   for (int b = 0; b < NZ; ++b)
      hRes[b] = new TH1F(Form("hRes_%d", b), "", 80, -8., 8.);

   TH1F *hResAll = new TH1F("hResAll", ";perp residual (mm);hits", 80, -8., 8.);
   TH1F *hResMean = new TH1F("hResMean", ";<perp residual per track> (mm);tracks", 80, -8., 8.);

   Long64_t n = tDigi->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      tDigi->GetEntry(i);
      if (patEvtArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patEvtArr->At(0);
      auto &tracks = pat->GetTrackCand();
      if (tracks.empty()) continue;
      const AtTrack *best = nullptr;
      size_t bestN = 0;
      for (const auto &t : tracks) {
         size_t nh = t.GetHitArray().size();
         if (nh > bestN) { bestN = nh; best = &t; }
      }
      if (!best || bestN < 10) continue;
      auto gc = best->GetGeoCenter();
      double cx = gc.first, cy = gc.second;
      double R = best->GetGeoRadius();
      if (R <= 0) continue;

      double sumRes = 0;
      int nh = 0;
      for (const auto &h : best->GetHitArray()) {
         const auto &p = h->GetPosition();
         double dx = p.X() - cx;
         double dy = p.Y() - cy;
         double res = std::sqrt(dx * dx + dy * dy) - R;
         hResAll->Fill(res);
         int b = -1;
         for (int k = 0; k < NZ; ++k)
            if (p.Z() >= zEdges[k] && p.Z() < zEdges[k + 1]) { b = k; break; }
         if (b >= 0) hRes[b]->Fill(res);
         sumRes += res;
         ++nh;
      }
      if (nh > 0) hResMean->Fill(sumRes / nh);
   }

   std::cout << "\n=== per-hit perpendicular residuals " << digiFile << " ===\n";
   std::cout << "z range (mm)    Nhits     mean (mm)   σ (mm)    (1/σ²=weight predicted)\n";
   std::cout << std::string(80, '-') << "\n";
   for (int b = 0; b < NZ; ++b) {
      double mean = hRes[b]->GetMean(), rms = hRes[b]->GetRMS();
      if (hRes[b]->GetEntries() > 30) {
         hRes[b]->Fit("gaus", "Q0", "", mean - 2.5 * rms, mean + 2.5 * rms);
         auto *f = hRes[b]->GetFunction("gaus");
         if (f) { mean = f->GetParameter(1); rms = f->GetParameter(2); }
      }
      printf("%4.0f-%-4.0f      %5.0f     %+5.2f      %5.2f\n",
             zEdges[b], zEdges[b + 1], hRes[b]->GetEntries(), mean, rms);
   }
   printf("ALL HITS        %5.0f     %+5.2f      %5.2f\n",
          hResAll->GetEntries(), hResAll->GetMean(), hResAll->GetRMS());
   printf("PER-TRACK MEAN  %5.0f     %+5.2f      %5.2f   <-- correlated noise indicator\n",
          hResMean->GetEntries(), hResMean->GetMean(), hResMean->GetRMS());
   std::cout << "\nInterpretation:\n";
   std::cout << " - If σ per z bin grows with √z → diffusion-limited (drift-aware should help).\n";
   std::cout << " - If σ is ~flat ≈ pad/√12 → granularity-limited.\n";
   std::cout << " - If PER-TRACK MEAN σ is comparable to ALL-HIT σ → correlated track-level\n";
   std::cout << "   noise dominates (whole track is shifted, weighting won't help).\n";
}
