/// @file write_gnn_digi.C
/// @brief Write GNN-separated hits into an AtPatternEvent ROOT file.
///
/// Creates output_gnn_digi.root where each AtTrack contains the raw hits
/// from one GNN-predicted track. The UKF's adaptive re-clustering
/// (ClusterizeSmooth3D) then clusters them — combining GNN's track
/// separation with Smooth3D's clustering quality.
///
/// Prerequisite: run gnn_hits_export.py first.
///
/// Run: root -b -q write_gnn_digi.C

#include <fstream>
#include <sstream>
#include <map>
#include <vector>

struct GNNHit {
   int hit_id, track_pred;
   double x, y, z, charge;
};

struct GNNTrackGeom {
   int track_pred, n_hits;
   double geo_radius, geo_theta, geo_phi, geo_cx, geo_cy;
};

std::map<int, std::vector<GNNHit>> LoadHits(const std::string &path)
{
   std::map<int, std::vector<GNNHit>> result;
   std::ifstream f(path);
   if (!f.is_open()) { std::cerr << "Cannot open " << path << std::endl; return result; }
   std::string line;
   std::getline(f, line); // header
   while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string tok;
      int eid;
      std::getline(ss, tok, ','); eid = std::stoi(tok);
      GNNHit h;
      std::getline(ss, tok, ','); h.hit_id = std::stoi(tok);
      std::getline(ss, tok, ','); h.track_pred = std::stoi(tok);
      std::getline(ss, tok, ','); h.x = std::stod(tok);
      std::getline(ss, tok, ','); h.y = std::stod(tok);
      std::getline(ss, tok, ','); h.z = std::stod(tok);
      std::getline(ss, tok, ','); h.charge = std::stod(tok);
      result[eid].push_back(h);
   }
   return result;
}

std::map<int, std::vector<GNNTrackGeom>> LoadTracks(const std::string &path)
{
   std::map<int, std::vector<GNNTrackGeom>> result;
   std::ifstream f(path);
   if (!f.is_open()) { std::cerr << "Cannot open " << path << std::endl; return result; }
   std::string line;
   std::getline(f, line);
   while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string tok;
      int eid;
      std::getline(ss, tok, ','); eid = std::stoi(tok);
      GNNTrackGeom tg;
      std::getline(ss, tok, ','); tg.track_pred = std::stoi(tok);
      std::getline(ss, tok, ','); tg.n_hits = std::stoi(tok);
      std::getline(ss, tok, ','); tg.geo_radius = std::stod(tok);
      std::getline(ss, tok, ','); tg.geo_theta = std::stod(tok);
      std::getline(ss, tok, ','); tg.geo_phi = std::stod(tok);
      std::getline(ss, tok, ','); tg.geo_cx = std::stod(tok);
      std::getline(ss, tok, ','); tg.geo_cy = std::stod(tok);
      result[eid].push_back(tg);
   }
   return result;
}

void write_gnn_digi(int nEvents = 5002)
{
   auto hits = LoadHits("data/gnn_training/gnn_separated_hits.csv");
   auto tracks = LoadTracks("data/gnn_training/gnn_separated_tracks.csv");
   std::cout << "Loaded GNN hits for " << hits.size() << " events" << std::endl;

   // Get event count from original digi file
   TFile *fIn = TFile::Open("data/output_digi.root");
   if (!fIn) { std::cerr << "Cannot open output_digi.root" << std::endl; return; }
   TTree *tIn = (TTree *)fIn->Get("cbmsim");
   int nEntries = std::min(nEvents, (int)tIn->GetEntries());

   // Create output with FairRoot metadata
   TFile *fOut = new TFile("data/output_gnn_digi.root", "RECREATE");

   // FairRoot requires these objects
   TFolder *cbmout = fOut->mkdir("cbmout") ? nullptr : nullptr;
   fOut->mkdir("cbmout");

   TList *branchList = new TList();
   branchList->SetName("BranchList");
   branchList->Add(new TObjString("AtPatternEvent"));

   TList *tbBranchList = new TList();
   tbBranchList->SetName("TimeBasedBranchList");

   fOut->cd();

   TTree *tOut = new TTree("cbmsim", "/cbmout");

   TClonesArray *patEvtArr = new TClonesArray("AtPatternEvent", 1);
   tOut->Branch("AtPatternEvent", &patEvtArr);

   int nTracks = 0;

   for (int iEv = 0; iEv < nEntries; iEv++) {
      patEvtArr->Clear("C");
      AtPatternEvent *patEvt = (AtPatternEvent *)patEvtArr->ConstructedAt(0);

      auto itH = hits.find(iEv);
      auto itT = tracks.find(iEv);

      if (itH != hits.end() && itT != tracks.end()) {
         // Group hits by track_pred
         std::map<int, std::vector<GNNHit *>> hitsByTrack;
         for (auto &h : itH->second)
            hitsByTrack[h.track_pred].push_back(&h);

         // Create one AtTrack per GNN track, with raw hits
         for (auto &tg : itT->second) {
            auto itHT = hitsByTrack.find(tg.track_pred);
            if (itHT == hitsByTrack.end() || (int)itHT->second.size() < 5)
               continue;

            // Select proton-like tracks using observables:
            //   theta_lab = 180 - theta_digi (geo_theta is in digi frame, degrees)
            //   Protons: theta_lab 20-90°, circle radius > 30 mm
            //   Beam: theta_lab < 20° or radius < 10 mm
            double theta_lab = 180.0 - tg.geo_theta;
            if (theta_lab < 20.0 || theta_lab > 90.0)
               continue;
            if (tg.geo_radius < 30.0)
               continue;

            AtTrack track;
            track.SetGeoRadius(tg.geo_radius);
            track.SetGeoTheta(tg.geo_theta * M_PI / 180.0);
            track.SetGeoPhi(tg.geo_phi * M_PI / 180.0);
            track.SetGeoCenter({tg.geo_cx, tg.geo_cy});

            // Add raw hits (NOT clusters) — ClusterizeSmooth3D will cluster them
            for (auto *h : itHT->second) {
               AtHit hit;
               hit.SetPosition({h->x, h->y, h->z});
               hit.SetCharge(h->charge);
               hit.SetTimeStamp(static_cast<int>(h->z / 10.0)); // approximate TB from Z
               track.AddHit(hit);
            }

            // Pre-cluster with ClusterizeSmooth3D so AtUKFDisplay can navigate
            AtTools::AtTrackTransformer transformer;
            transformer.ClusterizeSmooth3D(track, 20.0, 15.0);

            patEvt->AddTrack(track);
            nTracks++;
         }
      }

      tOut->Fill();

      if ((iEv + 1) % 1000 == 0)
         std::cout << "\r Writing " << iEv + 1 << "/" << nEntries << "..." << std::flush;
   }

   tOut->Write();
   branchList->Write("BranchList", TObject::kSingleKey);
   tbBranchList->Write("TimeBasedBranchList", TObject::kSingleKey);
   fOut->Close();
   fIn->Close();

   std::cout << "\n\nWrote " << nEntries << " events (" << nTracks << " tracks) to data/output_gnn_digi.root"
             << std::endl;
   std::cout << "\nNow run the UKF with adaptive re-clustering enabled:" << std::endl;
   std::cout << "  root -b -q run_ukf_gnn.C" << std::endl;
   std::cout << "Or visualize with the full display:" << std::endl;
   std::cout << "  root -l run_ukf_display_gnn.C" << std::endl;
}
