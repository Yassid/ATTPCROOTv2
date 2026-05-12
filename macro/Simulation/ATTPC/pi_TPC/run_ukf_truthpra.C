/// @file run_ukf_truthpra.C
/// @brief UKF on truth-seeded track candidates (bypass AtPRAtask).
///
/// For each event: take the digi-frame hits from AtEvent, build a single
/// AtTrack with all of them, and seed (GeoTheta, GeoPhi, GeoRadius,
/// GeoCenter) from MC truth. Run AtFitterUKF on this synthetic
/// AtPatternEvent and write the AtTrackingEvent to
/// data/output_ukf_truthpra.root. This isolates the *intrinsic* UKF
/// performance vs polar angle from the AtPRAtask candidate-finding
/// inefficiency we observed for backward (theta>100°) tracks.
///
/// Run: root -b -q 'run_ukf_truthpra.C(2000, +1)'

void run_ukf_truthpra(int nEvents = 2000, int pionSign = +1,
                      double momSigmaFrac = 0.1, int nIterations = 3,
                      const char *outSuffix = "", double eLossScale = 1.0)
{
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   // ---- physics constants ---------------------------------------------
   const double e_C = 1.602176634e-19;
   const double mass_pi_MeV = 139.57039;
   const double u_MeV = 931.49410372;
   const double mass_pi_amu = mass_pi_MeV / u_MeV;
   const double charge = (pionSign >= 0) ? +e_C : -e_C;
   const double Bz_T = 0.5;
   const double zPadPlane_mm = 1000.0;

   // ---- UKF setup (mirror run_ukf_only.C) -----------------------------
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(1.654e-3);
   eloss->SetProjectile(1, 1, mass_pi_amu);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(18, 40, 9));
   mat.push_back(std::make_tuple(6, 12, 1));
   mat.push_back(std::make_tuple(1, 1, 4));
   eloss->SetMaterial(mat);

   EventFit::AtFitterUKF ukf(charge, mass_pi_MeV, std::move(eloss));
   ukf.SetBField(ROOT::Math::XYZVector(0, 0, Bz_T));
   ukf.SetUKFParameters(1e-3, 2.0, 0.0);
   ukf.SetMeasurementSigma(1.0);
   ukf.SetMomentumSigmaFrac(momSigmaFrac);
   ukf.SetEnableEnergyStraggling(true);
   ukf.SetELossScaleFactor(eLossScale);
   ukf.SetMinClusters(5);
   ukf.SetNIterations(nIterations);
   ukf.SetZPadPlane(zPadPlane_mm);
   ukf.SetBackExtrapMaxPath(250.0);
   // Use circle-tangent seed (sort clusters by xy distance from beam axis)
   // instead of legacy chord seed. The chord seed assumes seed=track-start,
   // which is wrong for backward (theta_lab>90°) tracks where the
   // highest-Z_digi cluster is the track *end*, not the vertex side.
   ukf.SetUseClusterDirSeed(true);

   // ---- I/O -----------------------------------------------------------
   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *pts = new TClonesArray("AtMCPoint");
   TClonesArray *evArr = new TClonesArray("AtEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tSim->SetBranchAddress("AtTpcPoint", &pts);
   tDigi->SetBranchAddress("AtEventH", &evArr);

   TString outName = TString("data/output_ukf_truthpra") + outSuffix + ".root";
   TFile fOut(outName, "RECREATE");
   TTree tOut("cbmsim", "truth-PRA UKF output");
   TClonesArray teOutArr("AtTrackingEvent", 1);
   TClonesArray fmOutArr("AtFitMetadata", 1);
   tOut.Branch("AtTrackingEvent", &teOutArr);
   tOut.Branch("AtFitMetadata", &fmOutArr);

   const Long64_t n = std::min<Long64_t>(nEvents, std::min(tSim->GetEntries(), tDigi->GetEntries()));
   int nFit = 0, nBuilt = 0, nSkipNoMC = 0, nSkipNoHits = 0;

   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);

      teOutArr.Delete();
      fmOutArr.Delete();
      auto *te = (AtTrackingEvent *)teOutArr.ConstructedAt(0);
      auto *fm = (AtFitMetadata *)fmOutArr.ConstructedAt(0);

      // MC truth
      if (trks->GetEntries() == 0) { ++nSkipNoMC; tOut.Fill(); continue; }
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) { ++nSkipNoMC; tOut.Fill(); continue; }

      // Seed momentum/direction at the SEED CLUSTER LOCATION, which is the
      // highest-z_digi cluster. In digi frame z_digi = ZPadPlane - z_lab*10,
      // so highest z_digi == lowest z_lab. Find the MC point with min z_lab
      // (in the drift volume, primary track only) and read its momentum
      // there. This gives the correct R/phi at the seed: for forward tracks
      // this is the vertex (truth p_initial), for backward tracks this is
      // the far end (after energy loss, where the UKF actually starts).
      double zLabMin = 1e30;
      double pxS = 0, pyS = 0, pzS = 0;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetTrackID() != 0) continue;
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (p->GetZ() < zLabMin) {
            zLabMin = p->GetZ();
            pxS = p->GetPx();
            pyS = p->GetPy();
            pzS = p->GetPz();
         }
      }
      if (zLabMin > 1e29) { ++nSkipNoHits; tOut.Fill(); continue; }
      double pSeed_GeV = std::sqrt(pxS * pxS + pyS * pyS + pzS * pzS);
      double pSeed_MeV = pSeed_GeV * 1000.0;
      double thetaSeed = std::acos(pzS / pSeed_GeV); // rad, lab frame
      double phiSeed = std::atan2(pyS, pxS);
      double pTseed = pSeed_MeV * std::sin(thetaSeed);
      double R_mm = pTseed / (0.3 * Bz_T);

      // Helix center in xy at the seed point. For pi+ in +Bz, the center is
      // offset perpendicular to pT — the sign convention here matches the
      // initial-condition test that gave 2-8% σ/E for forward tracks.
      // We use the SEED MC point as the local position (vertex for forward,
      // far end for backward).
      double vx_mm = 0, vy_mm = 0;
      // Find the seed MC point's xy
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetTrackID() != 0) continue;
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (std::abs(p->GetZ() - zLabMin) < 1e-6) {
            vx_mm = p->GetX() * 10.0;
            vy_mm = p->GetY() * 10.0;
            break;
         }
      }
      double signQ = (pionSign >= 0) ? +1.0 : -1.0;
      double cx = vx_mm + signQ * R_mm * std::sin(phiSeed);
      double cy = vy_mm - signQ * R_mm * std::cos(phiSeed);

      // GeoTheta in digi frame: theta_lab = 180° - theta_digi · 180/π
      double geoTheta_digi = M_PI - thetaSeed;
      double geoPhi = phiSeed; // z-flip preserves phi

      // Hits from AtEvent (already in digi frame)
      if (evArr->GetEntries() == 0) { ++nSkipNoHits; tOut.Fill(); continue; }
      auto *ev = (AtEvent *)evArr->At(0);
      const auto &hits = ev->GetHits();
      if (hits.size() < 5) { ++nSkipNoHits; tOut.Fill(); continue; }

      AtPatternEvent pat;
      AtTrack track;
      for (const auto &h : hits)
         track.AddHit(*h);
      track.SetGeoTheta(geoTheta_digi);
      track.SetGeoPhi(geoPhi);
      track.SetGeoRadius(R_mm);
      track.SetGeoCenter(std::make_pair(cx, cy));
      pat.AddTrack(std::move(track));
      ++nBuilt;

      ukf.FitEvent(te, &pat, fm);
      auto &fitted = te->GetFittedTracks();
      if (!fitted.empty()) ++nFit;

      tOut.Fill();
      if ((i + 1) % 200 == 0)
         std::cout << "\r  Event " << (i + 1) << "/" << n
                   << "  built=" << nBuilt << "  fit=" << nFit << std::flush;
   }
   std::cout << "\n";

   fOut.cd();
   tOut.Write();
   fOut.Close();

   std::cout << "\n=== truth-PRA UKF run summary ===\n";
   std::cout << "Events processed:  " << n << "\n";
   std::cout << "Skipped (no MC):   " << nSkipNoMC << "\n";
   std::cout << "Skipped (no hits): " << nSkipNoHits << "\n";
   std::cout << "Tracks built:      " << nBuilt << "\n";
   std::cout << "Tracks fitted:     " << nFit << "\n";
   std::cout << "Wrote: " << outName << "\n";
}
