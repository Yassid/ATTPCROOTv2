/// Check what z range AtHits cover for forward vs backward tracks.
void inspect_hits()
{
   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *evArr = new TClonesArray("AtEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tDigi->SetBranchAddress("AtEventH", &evArr);
   tDigi->SetBranchAddress("AtPatternEvent", &patArr);

   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());

   int nForward = 0, nBackward = 0;
   for (Long64_t i = 0; i < n && (nForward < 3 || nBackward < 3); ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double pmc = std::sqrt(mc->GetPx()*mc->GetPx() + mc->GetPy()*mc->GetPy() + mc->GetPz()*mc->GetPz());
      double thMC = std::acos(mc->GetPz()/pmc) * 180. / M_PI;

      bool fwd = thMC < 70;
      bool bwd = thMC > 110 && thMC < 140;
      if (!fwd && !bwd) continue;
      if (fwd && nForward >= 3) continue;
      if (bwd && nBackward >= 3) continue;

      if (evArr->GetEntries() == 0) continue;
      auto *ev = (AtEvent *)evArr->At(0);
      const auto &hits = ev->GetHits();
      if (hits.size() < 5) continue;

      double zMin=1e9, zMax=-1e9, xMin=1e9, xMax=-1e9, yMin=1e9, yMax=-1e9;
      for (const auto &h : hits) {
         const auto &p = h->GetPosition();
         if (p.Z() < zMin) zMin = p.Z();
         if (p.Z() > zMax) zMax = p.Z();
         if (p.X() < xMin) xMin = p.X();
         if (p.X() > xMax) xMax = p.X();
         if (p.Y() < yMin) yMin = p.Y();
         if (p.Y() > yMax) yMax = p.Y();
      }

      // PRA candidate Geo*
      double R = NAN, gth = NAN, gph = NAN;
      int nPra = 0;
      if (patArr->GetEntries() > 0) {
         auto *pat = (AtPatternEvent *)patArr->At(0);
         nPra = pat->GetTrackCand().size();
         if (nPra > 0) {
            auto &t = pat->GetTrackCand()[0];
            R = t.GetGeoRadius();
            gth = t.GetGeoTheta() * 180./M_PI;
            gph = t.GetGeoPhi() * 180./M_PI;
         }
      }

      printf("Evt %lld  thMC=%.1f°  hits=%lu  z[%.1f,%.1f] x[%.1f,%.1f] y[%.1f,%.1f]  PRA: nCand=%d R=%.1f gth=%.1f° gph=%.1f°\n",
             i, thMC, hits.size(), zMin, zMax, xMin, xMax, yMin, yMax, nPra, R, gth, gph);
      if (fwd) ++nForward;
      if (bwd) ++nBackward;
   }
}
