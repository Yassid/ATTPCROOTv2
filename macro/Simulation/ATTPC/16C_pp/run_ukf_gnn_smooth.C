/// @file run_ukf_gnn_smooth.C
/// @brief Run UKF on GNN-separated tracks with ClusterizeSmooth3D clustering.
///
/// GNN for track separation + ClusterizeSmooth3D for clustering.
/// Reads output_gnn_digi.root (AtPatternEvent with raw hits per track),
/// runs FitEvent which triggers adaptive re-clustering internally.
///
/// Prerequisites:
///   1. Run gnn_hits_export.py → gnn_separated_hits.csv
///   2. Run write_gnn_digi.C → output_gnn_digi.root
///
/// Run: root -b -q run_ukf_gnn_smooth.C

void run_ukf_gnn_smooth(int nEvents = 5002)
{
   TStopwatch timer;
   timer.Start();

   // Open GNN digi file
   TFile *fIn = TFile::Open("data/output_gnn_digi.root");
   if (!fIn) { std::cerr << "Cannot open output_gnn_digi.root" << std::endl; return; }
   TTree *tIn = (TTree *)fIn->Get("cbmsim");

   TClonesArray *patEvtArr = nullptr;
   tIn->SetBranchAddress("AtPatternEvent", &patEvtArr);

   int nEntries = std::min(nEvents, (int)tIn->GetEntries());

   // Output
   TFile *fOut = new TFile("data/output_ukf_gnn_smooth.root", "RECREATE");
   TTree *tOut = new TTree("cbmsim", "UKF GNN+Smooth3D results");

   TClonesArray *trackingEvtArr = new TClonesArray("AtTrackingEvent", 1);
   tOut->Branch("AtTrackingEvent", &trackingEvtArr);

   // UKF fitter — same config as run_ukf_only.C
   double charge = 1.602176634e-19;
   double mass = 938.272;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
   eloss->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1));
   eloss->SetMaterial(mat);

   auto fitter = std::make_unique<EventFit::AtFitterUKF>(charge, mass, std::move(eloss));
   fitter->SetBField(ROOT::Math::XYZVector(0, 0, 2.85));
   fitter->SetUKFParameters(1e-3, 2.0, 0.0);
   fitter->SetMeasurementSigma(2.0);
   fitter->SetMomentumSigmaFrac(0.5);
   fitter->SetEnableEnergyStraggling(false);
   fitter->SetMinClusters(10);
   fitter->SetNIterations(1);
   fitter->SetZPadPlane(1000.0);
   // Adaptive re-clustering ON — re-cluster the raw hits with ClusterizeSmooth3D
   fitter->SetAdaptiveClustering(true);

   int nFitted = 0;

   for (int iEv = 0; iEv < nEntries; iEv++) {
      tIn->GetEntry(iEv);
      trackingEvtArr->Delete();

      if (!patEvtArr || patEvtArr->GetEntries() == 0) {
         trackingEvtArr->ConstructedAt(0);
         tOut->Fill();
         continue;
      }

      auto *patEvt = (AtPatternEvent *)patEvtArr->At(0);
      auto *trackingEvt = (AtTrackingEvent *)trackingEvtArr->ConstructedAt(0);

      fitter->FitEvent(trackingEvt, patEvt);

      auto &fittedTracks = trackingEvt->GetFittedTracks();
      nFitted += fittedTracks.size();

      if ((iEv + 1) % 500 == 0)
         std::cout << "\r Event " << iEv + 1 << "/" << nEntries
                   << " (" << nFitted << " tracks fitted)..." << std::flush;

      tOut->Fill();
   }

   tOut->Write();
   fOut->Close();
   fIn->Close();

   timer.Stop();

   std::cout << "\n\n=== UKF GNN+Smooth3D Results ===" << std::endl;
   std::cout << "Events: " << nEntries << std::endl;
   std::cout << "Tracks fitted: " << nFitted << std::endl;
   std::cout << "Output: data/output_ukf_gnn_smooth.root" << std::endl;
   std::cout << "Time: " << timer.RealTime() << " s" << std::endl;
}
