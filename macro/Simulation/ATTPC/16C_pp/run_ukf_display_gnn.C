/// @file run_ukf_display_gnn.C
/// @brief Launch AtUKFDisplay with GNN arc-walk clusters.
///
/// Same full interactive GUI as run_ukf_display.C (control panel, live
/// refitting, diagnostics) but reads GNN-clustered data instead of
/// RANSAC/Smooth3D clusters.
///
/// Prerequisites:
///   1. Run gnn_cluster_export.py → CSV files
///   2. Run write_gnn_digi.C → output_gnn_digi.root
///   3. Optionally run run_ukf_gnn.C → output_ukf_gnn.root
///
/// Run: root -l run_ukf_display_gnn.C
///      root -l 'run_ukf_display_gnn.C(3)' // start at event 3

void run_ukf_display_gnn(int startEvent = 1)
{
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   auto *display = AtUKFDisplay::GetInstance();

   // Load GNN-clustered events as "digi" input.
   // The display reads AtPatternEvent from this file — the clusters
   // are already GNN arc-walk clusters, so disable adaptive re-clustering
   // in the GUI (or set it via display options).
   display->LoadFiles(
      "data/output_gnn_digi.root",         // GNN-separated tracks with Smooth3D clusters
      "data/output_ukf_gnn_smooth.root",   // UKF fit results (AtTrackingEvent)
      "data/attpcsim.root"                 // MC truth (optional)
   );
   display->GotoEvent(startEvent);
}
