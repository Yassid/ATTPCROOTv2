/// @file extract_gnn_data.C
/// @brief Extract AT-TPC hit data with MC truth labels for GNN training.
///
/// Reads the digi output (AtEvent with hits + MC truth) and writes a flat
/// CSV file per event, plus a summary file. Each hit becomes one row with:
///   event_id, hit_id, x, y, z, charge, r, phi, track_id, A, Z, energy, eloss
///
/// Usage: root -b -q 'extract_gnn_data.C(10000)'

void extract_gnn_data(int nEvents = 10000)
{
   TFile *fd = TFile::Open("data/output_digi.root");
   if (!fd) { std::cerr << "Cannot open digi file" << std::endl; return; }
   TTree *td = (TTree *)fd->Get("cbmsim");

   TClonesArray *evtArr = nullptr;
   td->SetBranchAddress("AtEventH", &evtArr);

   int nEntries = std::min(nEvents, (int)td->GetEntries());

   // Create output directory
   gSystem->mkdir("data/gnn_training", kTRUE);

   // Single file with all hits
   std::ofstream fAll("data/gnn_training/all_hits.csv");
   fAll << "event_id,hit_id,x,y,z,charge,r,phi,track_id,A,Z,energy,eloss" << std::endl;

   // Summary file
   std::ofstream fSum("data/gnn_training/event_summary.csv");
   fSum << "event_id,n_hits,n_tracks,has_reaction" << std::endl;

   int nReaction = 0, nTotalHits = 0;

   for (int i = 0; i < nEntries; i++) {
      td->GetEntry(i);
      if (!evtArr || evtArr->GetEntries() == 0) continue;

      auto *evt = (AtEvent *)evtArr->At(0);
      int nHits = evt->GetNumHits();
      if (nHits == 0) continue;

      // Count tracks in this event
      std::set<int> trackIDs;
      bool hasMC = false;
      for (int h = 0; h < nHits; h++) {
         auto &hit = evt->GetHit(h);
         auto &mc = hit.GetMCSimPointArray();
         if (!mc.empty()) {
            hasMC = true;
            trackIDs.insert(mc[0].trackID);
         }
      }

      bool hasReaction = trackIDs.size() > 1; // beam + at least one reaction product
      if (hasReaction) nReaction++;

      fSum << i << "," << nHits << "," << trackIDs.size() << "," << (hasReaction ? 1 : 0) << std::endl;

      // Write hits
      for (int h = 0; h < nHits; h++) {
         auto &hit = evt->GetHit(h);
         auto pos = hit.GetPosition();
         double charge = hit.GetCharge();
         double r = std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y());
         double phi = std::atan2(pos.Y(), pos.X());

         int trackID = -1;
         int A = 0, Z = 0;
         double energy = 0, eloss = 0;
         auto &mc = hit.GetMCSimPointArray();
         if (!mc.empty()) {
            trackID = mc[0].trackID;
            A = mc[0].A;
            Z = mc[0].Z;
            energy = mc[0].energy;
            eloss = mc[0].eloss;
         }

         fAll << i << "," << h << ","
              << pos.X() << "," << pos.Y() << "," << pos.Z() << ","
              << charge << "," << r << "," << phi << ","
              << trackID << "," << A << "," << Z << ","
              << energy << "," << eloss << std::endl;
         nTotalHits++;
      }

      if ((i + 1) % 1000 == 0)
         std::cout << "\r Processed " << i + 1 << "/" << nEntries << " events..." << std::flush;
   }

   fAll.close();
   fSum.close();

   std::cout << "\n\n=== Extraction Summary ===" << std::endl;
   std::cout << "Events processed: " << nEntries << std::endl;
   std::cout << "Reaction events: " << nReaction << std::endl;
   std::cout << "Total hits: " << nTotalHits << std::endl;
   std::cout << "Files saved in data/gnn_training/" << std::endl;
   std::cout << "  all_hits.csv (" << nTotalHits << " rows)" << std::endl;
   std::cout << "  event_summary.csv" << std::endl;
}
