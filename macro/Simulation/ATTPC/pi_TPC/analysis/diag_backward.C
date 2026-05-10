/// @brief Diagnose why backward tracks (theta>90) have low fit yield.
/// Reads sim, digi, and ukf output; reports per theta bin: thrown, MC points,
/// pad hits (AtEvent), PRA candidates, fitted tracks.
void diag_backward()
{
   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   TFile fUKF("data/output_ukf_only.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *pts = new TClonesArray("AtMCPoint");
   TClonesArray *evArr = new TClonesArray("AtEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tSim->SetBranchAddress("AtTpcPoint", &pts);
   tDigi->SetBranchAddress("AtEventH", &evArr);
   tDigi->SetBranchAddress("AtPatternEvent", &patArr);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<int> nTh(NB, 0), nWithPts(NB, 0), nWithHits(NB, 0), nWithPra(NB, 0), nWithFit(NB, 0);
   std::vector<double> sumPts(NB, 0), sumHits(NB, 0), sumZmin(NB, 0), sumZmax(NB, 0);

   Long64_t n = std::min({tSim->GetEntries(), tDigi->GetEntries(), tUKF->GetEntries()});
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
      double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      int b = -1;
      for (int k = 0; k < NB; ++k)
         if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;

      int nPts = pts->GetEntries();
      ++nTh[b];
      sumPts[b] += nPts;
      double zMin = 1e9, zMax = -1e9;
      for (int j = 0; j < nPts; ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         double z = p->GetZ();
         if (z < zMin) zMin = z;
         if (z > zMax) zMax = z;
      }
      if (nPts > 0) { sumZmin[b] += zMin; sumZmax[b] += zMax; nWithPts[b]++; }

      int nHits = 0;
      if (evArr->GetEntries() > 0) {
         auto *ev = (AtEvent *)evArr->At(0);
         nHits = ev->GetNumHits();
      }
      sumHits[b] += nHits;
      if (nHits > 5) ++nWithHits[b];

      int nPra = 0;
      if (patArr->GetEntries() > 0) {
         auto *pat = (AtPatternEvent *)patArr->At(0);
         nPra = pat->GetTrackCand().size();
      }
      if (nPra > 0) ++nWithPra[b];

      int nFit = 0;
      if (teArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)teArr->At(0);
         nFit = te->GetFittedTracks().size();
      }
      if (nFit > 0) ++nWithFit[b];
   }

   std::cout << "\n=== Per-theta-bin attrition (sim->digi->PRA->UKF) ===\n";
   std::cout << "theta_MC      Nthr  <NMCpt> <zMin> <zMax>  hits>5  PRA   FIT  yield\n";
   std::cout << std::string(75, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      double avgPts = nTh[b] > 0 ? sumPts[b] / nTh[b] : 0;
      double avgZmin = nWithPts[b] > 0 ? sumZmin[b] / nWithPts[b] : 0;
      double avgZmax = nWithPts[b] > 0 ? sumZmax[b] / nWithPts[b] : 0;
      double yld = nTh[b] > 0 ? 100. * nWithFit[b] / nTh[b] : 0;
      char buf[256];
      snprintf(buf, sizeof(buf),
               "%4.0f-%4.0f    %4d   %5.0f   %5.1f  %5.1f    %4d  %4d  %4d  %5.1f%%",
               edges[b], edges[b + 1], nTh[b], avgPts, avgZmin, avgZmax,
               nWithHits[b], nWithPra[b], nWithFit[b], yld);
      std::cout << buf << "\n";
   }
}
