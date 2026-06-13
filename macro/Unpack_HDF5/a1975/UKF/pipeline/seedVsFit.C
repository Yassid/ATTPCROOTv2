/// @file seedVsFit.C
/// @brief For each fitted track, compare the PRA-derived SEED momentum (from the
/// source AtTrack's GeoRadius) against the genfit FITTED momentum, broken down by
/// theta region. Tells us whether collapsed/diverged KE come from a bad seed or a
/// diverging fit.  root -b -q 'pipeline/seedVsFit.C("/tmp/run_0106_genfitter_val.root")'

void seedVsFit(TString file = "/tmp/run_0106_genfitter_val.root", Double_t bField = 2.85, Double_t massAmu = 2.0135532)
{
   gSystem->Load("libAtReconstruction.so");
   TFile f(file);
   TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *trkEvt = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &trkEvt);
   const double mass = massAmu * 931.49410242;

   auto seedKE = [&](double radius_mm, double theta_rad) -> double {
      double sth = std::sin(theta_rad);
      if (std::abs(sth) < 1e-3) sth = 1e-3;
      double brho = bField * (radius_mm / 1000.0) / std::abs(sth);
      double p = 0.299792458 * brho * 1000.0; // MeV/c
      return std::sqrt(p * p + mass * mass) - mass;
   };

   const char *bandName[3] = {"forward th<45", "perp 70-110 ", "back th>135 "};
   long n[3] = {0, 0, 0}, seedBad[3] = {0, 0, 0}, fitCollapse[3] = {0, 0, 0}, fitDiverge[3] = {0, 0, 0};
   double seedSum[3] = {0, 0, 0}, fitSum[3] = {0, 0, 0};
   int printed = 0;

   Long64_t N = t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      auto *ev = (AtTrackingEvent *)trkEvt->At(0);
      if (!ev) continue;
      auto tracks = ev->GetTrackArray(); // source AtTracks (have GeoRadius)
      auto &fitted = ev->GetFittedTracks();
      for (auto &ftp : fitted) {
         auto *ft = ftp.get();
         if (!ft) continue;
         const auto &kin = ft->GetKinematics(0);
         double thdeg = kin.theta * TMath::RadToDeg();
         int b = -1;
         if (thdeg < 45) b = 0;
         else if (thdeg > 70 && thdeg < 110) b = 1;
         else if (thdeg > 135) b = 2;
         if (b < 0) continue;

         // find matching source track by ID -> its PRA radius
         double radius = -1;
         for (auto &tr : tracks)
            if (tr.GetTrackID() == ft->GetTrackID()) { radius = tr.GetGeoRadius(); break; }
         double sKE = (radius > 0) ? seedKE(radius, kin.theta) : -1;
         double fKE = kin.kineticEnergy;

         n[b]++;
         seedSum[b] += (sKE > 0 ? sKE : 0);
         fitSum[b] += fKE;
         if (sKE < 0 || sKE > 100) seedBad[b]++;
         if (fKE < 0.5) fitCollapse[b]++;
         if (fKE > 60) fitDiverge[b]++;

         if (b == 1 && printed < 25) { // sample the problematic perp blob
            printf("evt %5lld trk %2d  R_pra=%7.1fmm  seedKE=%8.2f  fitKE=%8.2f  th=%6.2f\n", i, ft->GetTrackID(),
                   radius, sKE, fKE, thdeg);
            printed++;
         }
      }
   }

   printf("\n=== seed vs fit by theta band ===\n");
   for (int b = 0; b < 3; ++b) {
      if (!n[b]) continue;
      printf("%s : N=%4ld  <seedKE>=%6.2f  <fitKE>=%6.2f  seedBad=%3ld(%.0f%%)  fitCollapse=%3ld(%.0f%%)  "
             "fitDiverge=%3ld(%.0f%%)\n",
             bandName[b], n[b], seedSum[b] / n[b], fitSum[b] / n[b], seedBad[b], 100.0 * seedBad[b] / n[b],
             fitCollapse[b], 100.0 * fitCollapse[b] / n[b], fitDiverge[b], 100.0 * fitDiverge[b] / n[b]);
   }
}
