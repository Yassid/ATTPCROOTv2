// Fit real tracks with MCMinimization::AtMCQMinimization and compare with the genfit fit of the
// same tracks. Written for run_0106 of a1975 (16C(p,p), H2 at 600 torr, 2.85 T), whose paths are
// the defaults below; point it at another run by passing them.
//
//   root -l -b -q 'MCQMinimizationOnData.C(20)'
//
// What it does with every single track event whose pattern is a circle:
//   - seeds the helix from that circle. The chord between the ends of the track is NOT the
//     initial direction once the track curves, so the dip angle comes from the swept arc against
//     the z extent and the azimuth is the tangent to the circle at the vertex,
//   - fits it and prints the result next to the genfit one, the distance between the end of the
//     simulated track and the last hit, the mean distance of each hit to the closest simulated
//     pad, and whether the simulated track curls the way the data does.
//
// NB: the sign of the magnetic field is not in the parameter file. On this run the tracks need
// bFieldSign = -1; with the wrong sign the simulated track curls away from the data and the mean
// distance goes from ~2 mm to ~35 mm.
//
// Reference result over the first 35 events of run_0106 (20 of them fitted by genfit as well):
// mean(this - genfit) = -1.0 deg in theta and +0.57 MeV in energy, with the simulated charge
// pattern sitting 0.1 to 8 mm from the hits.

void MCQMinimizationOnData(int maxTracks = 20, bool useRange = true, int bFieldSign = -1,
                           const char *patternFile = "/mnt/f/a1975/reco_pp/run_0106_multifit_reco.root",
                           const char *fitFile = "/mnt/f/a1975/reco/run_0106_genfitter_pphand.root",
                           const char *parFile = "ATTPC.a1954.par")
{
   TString repo = gSystem->Getenv("VMCWORKDIR");

   auto *run = new FairRunAna();
   auto *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open(repo + "/parameters/" + parFile, "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   rtdb->initContainers(0);

   auto map = std::make_shared<AtTpcMap>();
   map->ParseXMLMap(repo + "/scripts/Lookup20150611.xml");
   map->GeneratePadPlane();

   auto eLoss = std::make_shared<AtTools::AtELossCATIMA>(9.0e-5);
   eLoss->SetProjectile(1, 1, 938.272 / 931.494);
   std::vector<std::tuple<int, int, int>> mat;
   mat.emplace_back(1, 1, 1);
   eLoss->SetMaterial(mat);

   MCMinimization::AtMCQMinimization min;
   min.SetMap(map);
   min.SetParticle(1, 1);
   min.SetELossModel(eLoss);
   min.SetBackwardPropagation(false);
   min.SetUseRangeChi2(useRange);
   min.SetBFieldSign(bFieldSign);
   min.SetVerbose(false);

   TFile fPat(patternFile);
   auto *tPat = (TTree *)fPat.Get("cbmsim");
   TClonesArray *patArr = nullptr;
   tPat->SetBranchAddress("AtPatternEvent", &patArr);

   TFile fFit(fitFile);
   auto *tFit = (TTree *)fFit.Get("cbmsim");
   TClonesArray *fitArr = nullptr;
   tFit->SetBranchAddress("AtTrackingEvent", &fitArr);

   printf("\n entry nhits  seedTh   mcqTh  gfTh |   mcqE   gfE  |  simEnd-expEnd  meanDist  chi2      s   curl\n");
   printf("                 [deg]   [deg]  [deg]|  [MeV]  [MeV] |     [mm]          [mm]              \n");

   int done = 0;
   double sumdE = 0, sumdTh = 0;
   int nPaired = 0;
   for (Long64_t i = 0; i < tPat->GetEntries() && done < maxTracks; i++) {
      tPat->GetEntry(i);
      auto *pe = dynamic_cast<AtPatternEvent *>(patArr->At(0));
      if (!pe || pe->GetTrackCand().size() != 1)
         continue;
      auto &track = pe->GetTrackCand().at(0);
      const auto &hits = track.GetHitArray();
      if (hits.size() < 120)
         continue;
      const auto *circ = dynamic_cast<const AtPatterns::AtPatternCircle2D *>(track.GetPattern());
      if (!circ)
         continue;

      const AtHit *hv = hits.front().get(), *he = hits.front().get();
      for (const auto &h : hits) {
         if (h->GetPosition().Z() > hv->GetPosition().Z())
            hv = h.get();
         if (h->GetPosition().Z() < he->GetPosition().Z())
            he = h.get();
      }

      /* Seed the helix from the circle of the pattern recognition. The chord between the two ends
       * of the track is not the initial direction as soon as the track curves, so the swept
       * azimuth around the center of the circle is used instead:
       *    - the dip angle comes from the arc length against the z extent,
       *    - the azimuth is the tangent to the circle at the vertex. */
      const auto center = circ->GetCenter();
      const double R = circ->GetRadius();

      std::vector<std::pair<double, const AtHit *>> byZ; // hits from the vertex end downwards
      for (const auto &h : hits)
         byZ.emplace_back(-h->GetPosition().Z(), h.get());
      std::sort(byZ.begin(), byZ.end(), [](auto &a, auto &b) { return a.first < b.first; });

      // Average the hits in 20 slices along z, then unwrap the azimuth from slice to slice
      const int kSlices = 20;
      std::vector<double> sx(kSlices, 0), sy(kSlices, 0);
      std::vector<int> sn(kSlices, 0);
      const double zHi = hv->GetPosition().Z(), zLo = he->GetPosition().Z();
      for (const auto &e : byZ) {
         int s = int((zHi - e.second->GetPosition().Z()) / (zHi - zLo + 1e-9) * kSlices);
         s = std::max(0, std::min(kSlices - 1, s));
         sx[s] += e.second->GetPosition().X();
         sy[s] += e.second->GetPosition().Y();
         sn[s]++;
      }
      double sweep = 0, alphaPrev = 0, alphaFirst = 0;
      bool first = true;
      for (int s = 0; s < kSlices; s++) {
         if (sn[s] == 0)
            continue;
         const double a = std::atan2(sy[s] / sn[s] - center.Y(), sx[s] / sn[s] - center.X());
         if (first) {
            alphaFirst = alphaPrev = a;
            first = false;
            continue;
         }
         double d = a - alphaPrev;
         while (d > TMath::Pi())
            d -= TMath::TwoPi();
         while (d < -TMath::Pi())
            d += TMath::TwoPi();
         sweep += d;
         alphaPrev = a;
      }

      const double arc = std::abs(sweep) * R;            // path length in the pad plane [mm]
      const double dzTrack = std::abs(zHi - zLo);        // extent along the beam axis [mm]
      const double seedTheta = std::atan2(arc, dzTrack); // dip angle of the helix
      if (seedTheta * TMath::RadToDeg() < 30 || seedTheta * TMath::RadToDeg() > 89)
         continue;

      // Tangent to the circle at the vertex, pointing the way the track is swept
      const double sign = sweep > 0 ? 1. : -1.;
      const double tx = -sign * std::sin(alphaFirst);
      const double ty = sign * std::cos(alphaFirst);

      MCMinimization::AtMinimization::TrackSeed seed;
      seed.fVertex = ROOT::Math::XYZPoint(hv->GetPosition().X(), hv->GetPosition().Y(), hv->GetPosition().Z());
      seed.fVertexTB = hv->GetTimeStamp();
      seed.fTheta = seedTheta;
      seed.fPhi = TMath::Pi() - std::atan2(ty, tx);
      seed.fRadius = R;
      seed.fNumExpPoints = hits.size();

      TStopwatch watch;
      watch.Start();
      if (!min.Minimize(seed, track))
         continue;
      watch.Stop();
      const auto &res = min.GetFitPar();

      // Distance between the end of the simulated track and the last experimental hit
      const auto &simTrack = min.GetSimTrack();
      const double endDist = simTrack.empty() ? -1
                                              : std::sqrt(std::pow(simTrack.back().X() - he->GetPosition().X(), 2) +
                                                          std::pow(simTrack.back().Y() - he->GetPosition().Y(), 2) +
                                                          std::pow(simTrack.back().Z() - he->GetPosition().Z(), 2));

      // Mean distance of each hit to the closest pad of the simulated track
      const auto &simPads = min.GetSimPads();
      double meanDist = 0;
      for (const auto &h : hits) {
         double best = 1e9;
         for (const auto &s : simPads)
            best = std::min(best, std::hypot(h->GetPosition().X() - s.X(), h->GetPosition().Y() - s.Y()));
         meanDist += best;
      }
      meanDist /= hits.size();

      // Handedness: signed area swept by the simulated track, to compare with the data's sweep
      double simSweep = 0;
      for (size_t k = 2; k < simTrack.size(); k += 10) {
         const double ax = simTrack[k - 1].X() - simTrack[k - 2].X();
         const double ay = simTrack[k - 1].Y() - simTrack[k - 2].Y();
         const double bx = simTrack[k].X() - simTrack[k - 1].X();
         const double by = simTrack[k].Y() - simTrack[k - 1].Y();
         simSweep += ax * by - ay * bx;
      }

      // The genfit fit of the same event
      double gfE = -1, gfTh = -1;
      tFit->GetEntry(i);
      auto *te = dynamic_cast<AtTrackingEvent *>(fitArr->At(0));
      if (te && te->GetFittedTracks().size() == 1) {
         const auto &k = te->GetFittedTracks().at(0)->GetKinematics();
         gfE = k.kineticEnergy;
         gfTh = k.theta * TMath::RadToDeg();
      }

      printf("%6lld %5zu %7.1f %7.1f %6.1f | %6.3f %6.3f | %8.1f %10.1f %10.2e %6.1f   %s\n", i, hits.size(),
             seedTheta * TMath::RadToDeg(), res.fTheta * TMath::RadToDeg(), gfTh, res.fEnergy, gfE, endDist, meanDist,
             res.fChi2, watch.RealTime(), (sweep * simSweep > 0 ? "same" : "OPPOSITE"));
      if (gfE > 0) {
         sumdE += res.fEnergy - gfE;
         sumdTh += res.fTheta * TMath::RadToDeg() - gfTh;
         nPaired++;
      }
      done++;
   }

   if (nPaired)
      printf("\npaired with genfit: %d tracks, mean(mcq - genfit): energy %+.3f MeV, theta %+.1f deg\n", nPaired,
             sumdE / nPaired, sumdTh / nPaired);
}
