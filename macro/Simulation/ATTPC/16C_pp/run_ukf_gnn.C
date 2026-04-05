/// @file run_ukf_gnn.C
/// @brief Run UKF fitter with GNN arc-walk clusters instead of ClusterizeSmooth3D.
///
/// Reads output_digi.root (AtPatternEvent from RANSAC) and gnn_clusters.csv
/// (pre-computed by gnn_tracking/gnn_cluster_export.py), replaces cluster arrays,
/// and fits with the UKF.
///
/// Prerequisite:
///   cd gnn_tracking && python gnn_cluster_export.py && cd ..
///
/// Run: root -b -q run_ukf_gnn.C

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <memory>

// ─── CSV data structures ───────────────────────────────────────────

struct GNNCluster {
   int cluster_id;
   int track_pred;
   int n_hits;
   double x, y, z;
   double charge;
   double cov_xx, cov_yy, cov_zz;
};

struct GNNTrackGeom {
   int track_pred;
   int n_hits;
   int n_clusters;
   double geo_radius;
   double geo_theta;  // degrees
   double geo_phi;    // degrees
   double geo_cx, geo_cy;
};

// ─── CSV readers ───────────────────────────────────────────────────

std::map<int, std::vector<GNNCluster>> LoadGNNClusters(const std::string &path)
{
   std::map<int, std::vector<GNNCluster>> result;
   std::ifstream f(path);
   if (!f.is_open()) {
      std::cerr << "ERROR: Cannot open " << path << std::endl;
      return result;
   }

   std::string line;
   std::getline(f, line); // header

   while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string token;
      int event_id;

      std::getline(ss, token, ','); event_id = std::stoi(token);

      GNNCluster cl;
      std::getline(ss, token, ','); cl.cluster_id = std::stoi(token);
      std::getline(ss, token, ','); cl.track_pred = std::stoi(token);
      std::getline(ss, token, ','); cl.n_hits = std::stoi(token);
      std::getline(ss, token, ','); cl.x = std::stod(token);
      std::getline(ss, token, ','); cl.y = std::stod(token);
      std::getline(ss, token, ','); cl.z = std::stod(token);
      std::getline(ss, token, ','); cl.charge = std::stod(token);
      std::getline(ss, token, ','); cl.cov_xx = std::stod(token);
      std::getline(ss, token, ','); cl.cov_yy = std::stod(token);
      std::getline(ss, token, ','); cl.cov_zz = std::stod(token);

      result[event_id].push_back(cl);
   }

   std::cout << "Loaded GNN clusters for " << result.size() << " events" << std::endl;
   return result;
}

std::map<int, std::vector<GNNTrackGeom>> LoadGNNTracks(const std::string &path)
{
   std::map<int, std::vector<GNNTrackGeom>> result;
   std::ifstream f(path);
   if (!f.is_open()) {
      std::cerr << "ERROR: Cannot open " << path << std::endl;
      return result;
   }

   std::string line;
   std::getline(f, line); // header

   while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string token;
      int event_id;

      std::getline(ss, token, ','); event_id = std::stoi(token);

      GNNTrackGeom tg;
      std::getline(ss, token, ','); tg.track_pred = std::stoi(token);
      std::getline(ss, token, ','); tg.n_hits = std::stoi(token);
      std::getline(ss, token, ','); tg.n_clusters = std::stoi(token);
      std::getline(ss, token, ','); tg.geo_radius = std::stod(token);
      std::getline(ss, token, ','); tg.geo_theta = std::stod(token);
      std::getline(ss, token, ','); tg.geo_phi = std::stod(token);
      std::getline(ss, token, ','); tg.geo_cx = std::stod(token);
      std::getline(ss, token, ','); tg.geo_cy = std::stod(token);

      result[event_id].push_back(tg);
   }

   std::cout << "Loaded GNN track geometry for " << result.size() << " events" << std::endl;
   return result;
}

// ─── Main ──────────────────────────────────────────────────────────

void run_ukf_gnn(int nEvents = 10000)
{
   TStopwatch timer;
   timer.Start();

   // ── 1. Load GNN data ──
   auto gnnClusters = LoadGNNClusters("data/gnn_training/gnn_clusters.csv");
   auto gnnTracks = LoadGNNTracks("data/gnn_training/gnn_tracks.csv");

   if (gnnClusters.empty()) {
      std::cerr << "No GNN clusters loaded. Run gnn_cluster_export.py first." << std::endl;
      return;
   }

   // ── 2. Open input file ──
   TFile *fIn = TFile::Open("data/output_digi.root");
   if (!fIn) { std::cerr << "Cannot open data/output_digi.root" << std::endl; return; }
   TTree *tIn = (TTree *)fIn->Get("cbmsim");
   int nEntries = std::min(nEvents, (int)tIn->GetEntries());

   // We don't actually need AtPatternEvent from the file —
   // we build our own from GNN data. But we need the TTree
   // entry count to match event IDs.

   // ── 3. Create output file ──
   TFile *fOut = new TFile("data/output_ukf_gnn.root", "RECREATE");
   TTree *tOut = new TTree("convergence", "UKF GNN results");

   // Output branches — flat for easy analysis
   double fit_KE, fit_theta, fit_phi, fit_vtxX, fit_vtxY, fit_vtxZ;
   double fit_length;
   int fit_nClusters, fit_trackPred, out_event_id;
   int fit_status; // 0 = success, 1 = failed

   tOut->Branch("event_id", &out_event_id);
   tOut->Branch("track_pred", &fit_trackPred);
   tOut->Branch("KE", &fit_KE);
   tOut->Branch("theta", &fit_theta);
   tOut->Branch("phi", &fit_phi);
   tOut->Branch("vtxX", &fit_vtxX);
   tOut->Branch("vtxY", &fit_vtxY);
   tOut->Branch("vtxZ", &fit_vtxZ);
   tOut->Branch("length", &fit_length);
   tOut->Branch("nClusters", &fit_nClusters);
   tOut->Branch("status", &fit_status);

   // ── 4. Create UKF fitter ──
   double charge_C = 1.602176634e-19;
   double mass = 938.272;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1));
   eloss->SetMaterial(mat);

   auto fitter = std::make_unique<EventFit::AtFitterUKF>(charge_C, mass, std::move(eloss));
   fitter->SetBField(ROOT::Math::XYZVector(0, 0, 2.85));
   fitter->SetUKFParameters(1e-3, 2.0, 0.0);
   fitter->SetMeasurementSigma(2.0);
   fitter->SetMomentumSigmaFrac(0.5);
   fitter->SetEnableEnergyStraggling(false);
   fitter->SetMinClusters(10);
   fitter->SetNIterations(1);
   fitter->SetZPadPlane(1000.0);
   fitter->SetAdaptiveClustering(false); // GNN clusters are pre-computed

   int nFitted = 0, nSkipped = 0, nNoGNN = 0;

   // ── 5. Event loop ──
   for (int iEv = 0; iEv < nEntries; iEv++) {

      // Check if GNN data exists for this event
      auto itCl = gnnClusters.find(iEv);
      auto itTr = gnnTracks.find(iEv);
      if (itCl == gnnClusters.end() || itTr == gnnTracks.end()) {
         nNoGNN++;
         continue;
      }

      auto &evClusters = itCl->second;
      auto &evTracks = itTr->second;

      // Group GNN clusters by track_pred
      std::map<int, std::vector<GNNCluster *>> clustersByTrack;
      for (auto &cl : evClusters)
         clustersByTrack[cl.track_pred].push_back(&cl);

      // Build an AtPatternEvent with GNN-predicted tracks
      AtPatternEvent patEvt;

      for (auto &tg : evTracks) {
         int tp = tg.track_pred;
         auto itClT = clustersByTrack.find(tp);
         if (itClT == clustersByTrack.end() || (int)itClT->second.size() < 10)
            continue;

         AtTrack track;
         track.SetGeoRadius(tg.geo_radius);
         track.SetGeoTheta(tg.geo_theta * M_PI / 180.0);
         track.SetGeoPhi(tg.geo_phi * M_PI / 180.0);
         track.SetGeoCenter({tg.geo_cx, tg.geo_cy});

         // Add GNN clusters (already ordered by arc-walk: vertex end first)
         int clID = 0;
         for (auto *cl : itClT->second) {
            auto hitCluster = std::make_shared<AtHitCluster>();
            hitCluster->SetClusterID(clID++);
            hitCluster->SetPosition({cl->x, cl->y, cl->z});
            hitCluster->SetCharge(cl->charge);
            hitCluster->SetTimeStamp(0);

            TMatrixDSym cov(3);
            cov(0, 0) = cl->cov_xx;
            cov(1, 1) = cl->cov_yy;
            cov(2, 2) = cl->cov_zz;
            hitCluster->SetCovMatrix(cov);

            track.AddClusterHit(hitCluster);
         }

         patEvt.AddTrack(track);
      }

      if (patEvt.GetTrackCand().empty())
         continue;

      // Fit all tracks in this event via FitEvent
      AtTrackingEvent trackingEvent;
      fitter->FitEvent(&trackingEvent, &patEvt);

      // Extract results
      auto &fittedTracks = trackingEvent.GetFittedTracks();

      // We need to match fitted tracks back to GNN track_pred IDs.
      // FitEvent processes tracks in order, so fitted track i corresponds
      // to the i-th track candidate that was successfully fitted.
      // For simplicity, iterate over fitted tracks and record results.
      int trackIdx = 0;
      auto &trackCands = patEvt.GetTrackCand();

      for (auto &ft : fittedTracks) {
         out_event_id = iEv;
         // Try to identify which GNN track this is
         // Use the first candidate track index as approximation
         fit_trackPred = trackIdx < (int)evTracks.size() ? evTracks[trackIdx].track_pred : -1;

         auto &kin = ft->GetKinematics();
         fit_KE = kin.kineticEnergy;
         fit_theta = kin.theta;
         fit_phi = kin.phi;
         auto vtx = ft->GetVertex();
         fit_vtxX = vtx.X();
         fit_vtxY = vtx.Y();
         fit_vtxZ = vtx.Z();
         auto &props = ft->GetTrackPropertiesStruct();
         fit_length = props.trackLength;
         fit_nClusters = props.trackPoints;
         fit_status = 0;
         nFitted++;

         tOut->Fill();
         trackIdx++;
      }

      // Record skipped tracks
      int nCands = patEvt.GetTrackCand().size();
      nSkipped += nCands - (int)fittedTracks.size();

      if ((iEv + 1) % 1000 == 0)
         std::cout << "\r Processed " << iEv + 1 << "/" << nEntries << " events..." << std::flush;
   }

   tOut->Write();
   fOut->Close();
   fIn->Close();

   timer.Stop();

   std::cout << "\n\n=== UKF GNN Results ===" << std::endl;
   std::cout << "Events processed: " << nEntries << std::endl;
   std::cout << "Tracks fitted: " << nFitted << std::endl;
   std::cout << "Tracks skipped (fit failed): " << nSkipped << std::endl;
   std::cout << "Events without GNN data: " << nNoGNN << std::endl;
   std::cout << "Output: data/output_ukf_gnn.root" << std::endl;
   std::cout << "Time: " << timer.RealTime() << " s" << std::endl;
}
