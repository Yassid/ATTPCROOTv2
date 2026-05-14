/// @file diag_z_bias.C
/// @brief Diagnose the z-axis offset seen in display_event3d.
///
/// For one or more events, compare:
///   (a) MC point z values vs MC track Pz / propagation expectation
///   (b) Pad hit z values reconstructed via PSA
///   (c) UKF smoothed z values
///
/// Specifically check whether:
///   - The Geant4 transport bends in z (Bz ≠ 0 along expected direction?)
///   - The clusterize/pulse uses the same drift conv. as PSA's CalculateZGeo
///   - The MCPP target pad mapping shifts z when going through ToA(t)→z

void diag_z_bias(int P_MeV = 800, int nEvts = 5, const char *runDir = "data")
{
   TString simF = Form("%s/HYDRAsim_p%d.root", runDir, P_MeV);
   TString digiF = Form("%s/output_digi_p%d.root", runDir, P_MeV);
   TString ukfF  = Form("%s/output_ukf_HYDRA_p%d.root", runDir, P_MeV);

   TFile fS(simF), fD(digiF), fU(ukfF);
   auto *tS = (TTree *)fS.Get("cbmsim");
   auto *tD = (TTree *)fD.Get("cbmsim");
   auto *tU = (TTree *)fU.Get("cbmsim");

   auto *trks = new TClonesArray("AtMCTrack");
   auto *pts  = new TClonesArray("AtMCPoint");
   auto *ev   = new TClonesArray("AtEvent");
   auto *te   = new TClonesArray("AtTrackingEvent");
   tS->SetBranchAddress("MCTrack", &trks);
   tS->SetBranchAddress("AtTpcPoint", &pts);
   tD->SetBranchAddress("AtEventH", &ev);
   tU->SetBranchAddress("AtTrackingEvent", &te);

   for (int ie = 0; ie < nEvts; ++ie) {
      tS->GetEntry(ie); tD->GetEntry(ie); tU->GetEntry(ie);
      printf("\n========== Event %d ==========\n", ie);
      if (trks->GetEntries() == 0) { printf("no MC track\n"); continue; }

      auto *mc = (AtMCTrack *)trks->At(0);
      double vx = mc->GetStartX()*10, vy = mc->GetStartY()*10, vz = mc->GetStartZ()*10;
      double Px = mc->GetPx(), Py = mc->GetPy(), Pz = mc->GetPz();
      double p = std::sqrt(Px*Px + Py*Py + Pz*Pz);
      printf("MC vertex (mm): (%.2f, %.2f, %.2f)\n", vx, vy, vz);
      printf("MC momentum (GeV/c): (%.4f, %.4f, %.4f) |p|=%.4f, Pz/|p|=%.4e\n",
             Px, Py, Pz, p, Pz/p);

      // First, last MC point in drift_volume + range of z
      double zMinMC = 1e9, zMaxMC = -1e9;
      double xFirst = 1e9, zFirst = 0, yFirst = 0;
      int nMC = 0;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *pp = (AtMCPoint *)pts->At(j);
         if (pp->GetVolName() != TString("drift_volume")) continue;
         if (pp->GetTrackID() != 0) continue;
         double xm = pp->GetX()*10, ym = pp->GetY()*10, zm = pp->GetZ()*10;
         nMC++;
         zMinMC = std::min(zMinMC, zm);
         zMaxMC = std::max(zMaxMC, zm);
         if (xm < xFirst) {
            xFirst = xm; yFirst = ym; zFirst = zm;
         }
      }
      printf("MC pts in drift: %d, z range [%.3f, %.3f] mm (Δz = %.3f mm)\n",
             nMC, zMinMC, zMaxMC, zMaxMC - zMinMC);
      printf("MC first (in chamber, min x):  x=%.2f y=%.2f z=%.2f mm\n",
             xFirst, yFirst, zFirst);

      // Pad hits
      double zMinPad = 1e9, zMaxPad = -1e9, zSumPad = 0;
      int nPad = 0;
      double zPad_at_xFirst = 0; double distMin = 1e9;
      if (ev->GetEntries() > 0) {
         for (auto &h : ((AtEvent *)ev->At(0))->GetHits()) {
            auto pos = h->GetPosition();
            zMinPad = std::min(zMinPad, pos.Z());
            zMaxPad = std::max(zMaxPad, pos.Z());
            zSumPad += pos.Z(); nPad++;
            double d2 = (pos.X()-xFirst)*(pos.X()-xFirst) + (pos.Y()-yFirst)*(pos.Y()-yFirst);
            if (d2 < distMin) { distMin = d2; zPad_at_xFirst = pos.Z(); }
         }
      }
      printf("Pad hits: %d, z range [%.3f, %.3f] mm, mean=%.3f mm\n",
             nPad, zMinPad, zMaxPad, nPad?zSumPad/nPad:0);
      printf("Pad hit nearest MC-first (xy): z=%.3f mm   (Δz vs MC first = %.3f mm)\n",
             zPad_at_xFirst, zPad_at_xFirst - zFirst);

      // UKF
      if (te->GetEntries() > 0) {
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr; double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (best) {
            auto v = best->GetVertex();
            printf("UKF vertex (POCA): x=%.2f y=%.2f z=%.2f mm\n", v.X(), v.Y(), v.Z());
            double zMinU = 1e9, zMaxU = -1e9, zUatFirst = 0; double dMinU = 1e9;
            for (auto &pp : best->GetSmoothedPositions()) {
               zMinU = std::min(zMinU, pp.Z());
               zMaxU = std::max(zMaxU, pp.Z());
               double d2 = (pp.X()-xFirst)*(pp.X()-xFirst) + (pp.Y()-yFirst)*(pp.Y()-yFirst);
               if (d2 < dMinU) { dMinU = d2; zUatFirst = pp.Z(); }
            }
            printf("UKF smoothed z range [%.3f, %.3f] mm\n", zMinU, zMaxU);
            printf("UKF smoothed nearest MC-first (xy): z=%.3f mm  (Δz vs MC first = %.3f mm)\n",
                   zUatFirst, zUatFirst - zFirst);
         }
      }

      // Quick MC consistency: vertex at (-400, 44, 147) mm — what does Pz/p say?
      double vxm = vx, vzm = vz;
      printf("Vertex Pz = %.4f GeV (~ %.1f%% of |p|)  →  expected drift slope dz/dx ≈ %.3f\n",
             Pz, 100*Pz/p, Pz/std::sqrt(Px*Px+Py*Py));
   }
}
