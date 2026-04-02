/// @file run_ukf_digi.C
/// @brief End-to-end UKF validation with digitized data.
///
/// Reads the simulation (attpcsim.root) and digitized (output_digi.root) files,
/// extracts the proton track from reaction events, runs the UKF fitter on the
/// hit clusters, and compares the result against MC truth.
///
/// Prerequisites:
///   1. Run C16_pp_sim.C to generate attpcsim.root
///   2. Run run_digi_attpc.C to generate output_digi.root
///
/// Run: root -b -q run_ukf_digi.C
///      root -b -q 'run_ukf_digi.C(5)'          // process only 5 events
///      root -b -q 'run_ukf_digi.C(-1, 0.95)'   // apply energy loss scaling
///      root -b -q 'run_ukf_digi.C(-1, 1.0, 3.0, true)' // per-cluster covariance

using ROOT::Math::Polar3DVector;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

const double mass_p = 938.272;           // MeV/c^2
const double charge_p = 1.602176634e-19; // C
const double gasDensity = 3.553e-5;      // g/cm^3 for 300 torr H2
const double ZPadPlane = 1000.0;         // mm — pad plane Z position

/// Convert digi coordinates to MC (lab) frame.
/// X, Y are the same; Z_lab = ZPadPlane - Z_digi
XYZPoint DigiToLab(const XYZPoint &p)
{
   return {p.X(), p.Y(), ZPadPlane - p.Z()};
}

void run_ukf_digi(int maxEvents = -1, double eLossScale = 1.0, double minSpacing = 3.0,
                  bool usePerClusterCov = false)
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString inOutDir = "./data/";
   TString mcFile = inOutDir + "attpcsim.root";
   TString digiFile = inOutDir + "output_digi.root";

   // --- Open files ---
   auto fMC = TFile::Open(mcFile);
   auto fDigi = TFile::Open(digiFile);
   if (!fMC || !fDigi) {
      std::cerr << "Cannot open input files. Run C16_pp_sim.C and run_digi_attpc.C first." << std::endl;
      return;
   }

   auto treeMC = (TTree *)fMC->Get("cbmsim");
   auto treeDigi = (TTree *)fDigi->Get("cbmsim");

   TClonesArray *mcTracks = nullptr;
   TClonesArray *mcPoints = nullptr;
   treeMC->SetBranchAddress("MCTrack", &mcTracks);
   treeMC->SetBranchAddress("AtTpcPoint", &mcPoints);

   TClonesArray *patEvtArr = nullptr;
   treeDigi->SetBranchAddress("AtPatternEvent", &patEvtArr);

   int nEvents = treeDigi->GetEntries();
   if (maxEvents > 0 && maxEvents < nEvents)
      nEvents = maxEvents;

   // --- Create UKF ---
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   eloss->SetProjectile(1, 1, 1); // proton
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 1}); // H target
   eloss->SetMaterial(mat);

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(eloss));
   propagator.SetBField({0, 0, 2.85});
   auto stepper = std::make_unique<AtTools::AtRK4AdaptiveStepper>();
   kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));
   ukf.setParameters(1e-3, 2.0, 0.0);
   ukf.fEnableEnStraggling = true;
   ukf.fELossScaleFactor = eLossScale; // Scale CATIMA dE/dx (1.0 = no correction)

   if (std::abs(eLossScale - 1.0) > 1e-6)
      std::cout << "Energy loss scaling factor: " << eLossScale << std::endl;

   // --- Summary histograms ---
   TH1F *hMomErr = new TH1F("hMomErr", "Momentum Error;(p_{reco} - p_{true}) / p_{true} [%];Events", 50, -20, 20);
   TH1F *hThetaErr =
      new TH1F("hThetaErr", "#theta Error;#theta_{reco} - #theta_{true} [deg];Events", 50, -10, 10);
   TH1F *hResid = new TH1F("hResid", "Mean Cluster Residual;Residual [mm];Events", 50, 0, 10);
   TH1F *hNClusters = new TH1F("hNClusters", "Clusters per Track;N clusters;Events", 50, 0, 100);

   int nFitted = 0;
   int nFailed = 0;

   // --- Event loop ---
   for (int iEvt = 0; iEvt < nEvents; iEvt++) {
      treeMC->GetEntry(iEvt);
      treeDigi->GetEntry(iEvt);

      // Skip events without pattern recognition output
      if (patEvtArr->GetEntries() == 0)
         continue;
      auto *patEvt = (AtPatternEvent *)patEvtArr->At(0);
      auto &tracks = patEvt->GetTrackCand();
      if (tracks.empty())
         continue;

      // Find the proton MC track (PDG 2212) and collect MC points along it
      int protonTrackID = -1;
      XYZVector trueMomVertex; // momentum at vertex
      XYZPoint trueVertex;
      for (int j = 0; j < mcTracks->GetEntries(); j++) {
         auto *t = (AtMCTrack *)mcTracks->At(j);
         if (t->GetPdgCode() == 2212) {
            protonTrackID = j;
            trueMomVertex = XYZVector(t->GetPx() * 1e3, t->GetPy() * 1e3, t->GetPz() * 1e3); // GeV → MeV
            trueVertex = XYZPoint(t->GetStartX() * 10, t->GetStartY() * 10, t->GetStartZ() * 10); // cm → mm
            break;
         }
      }
      if (protonTrackID < 0)
         continue; // No proton in this event (beam-only event)

      // Collect proton MC points (position + momentum in lab frame, mm and MeV/c)
      struct MCPointInfo {
         XYZPoint pos;
         XYZVector mom;
      };
      std::vector<MCPointInfo> mcProtonPoints;
      for (int j = 0; j < mcPoints->GetEntries(); j++) {
         auto *p = (AtMCPoint *)mcPoints->At(j);
         if (p->GetTrackID() != protonTrackID)
            continue;
         mcProtonPoints.push_back(
            {{p->GetX() * 10, p->GetY() * 10, p->GetZ() * 10}, {p->GetPx() * 1e3, p->GetPy() * 1e3, p->GetPz() * 1e3}});
      }

      // Take the largest track (most hits = likely the proton)
      int bestTrack = 0;
      for (size_t t = 1; t < tracks.size(); t++) {
         if (tracks[t].GetHitArray().size() > tracks[bestTrack].GetHitArray().size())
            bestTrack = t;
      }

      auto &track = tracks[bestTrack];
      auto *clusters = track.GetHitClusterArray();
      if (clusters->size() < 5)
         continue;

      hNClusters->Fill(clusters->size());

      // --- Order clusters along the track arc and filter ---
      // Sort by Z first, then walk along the track greedily picking the
      // nearest-neighbor cluster that hasn't been visited, starting from
      // the cluster closest to the MC vertex.
      track.SortClusterHitArrayZ();

      // Convert all cluster positions to lab frame
      std::vector<XYZPoint> clusterLab(clusters->size());
      for (size_t i = 0; i < clusters->size(); i++)
         clusterLab[i] = DigiToLab(clusters->at(i).GetPosition());

      // Find the cluster closest to the MC vertex
      double minDist = 1e9;
      int seedIdx = 0;
      for (size_t i = 0; i < clusterLab.size(); i++) {
         double d = (clusterLab[i] - trueVertex).R();
         if (d < minDist) {
            minDist = d;
            seedIdx = i;
         }
      }

      // Nearest-neighbor walk to order along the arc (in lab frame)
      std::vector<int> order;
      std::vector<bool> used(clusterLab.size(), false);
      order.push_back(seedIdx);
      used[seedIdx] = true;
      for (size_t step = 1; step < clusterLab.size(); step++) {
         XYZPoint current = clusterLab[order.back()];
         double bestDist = 1e9;
         int bestIdx = -1;
         for (size_t j = 0; j < clusterLab.size(); j++) {
            if (used[j])
               continue;
            double d = (clusterLab[j] - current).R();
            if (d < bestDist) {
               bestDist = d;
               bestIdx = j;
            }
         }
         if (bestIdx < 0)
            break;
         order.push_back(bestIdx);
         used[bestIdx] = true;
      }

      // Build ordered cluster list in lab frame, skipping clusters too close.
      // Keep track of original cluster indices for per-cluster covariance.
      double minClusterSpacing = minSpacing;
      std::vector<XYZPoint> orderedClusters;
      std::vector<int> orderedIndices; // index into original clusters array
      orderedClusters.push_back(clusterLab[order[0]]);
      orderedIndices.push_back(order[0]);
      for (size_t i = 1; i < order.size(); i++) {
         XYZPoint pos = clusterLab[order[i]];
         if ((pos - orderedClusters.back()).R() >= minClusterSpacing) {
            orderedClusters.push_back(pos);
            orderedIndices.push_back(order[i]);
         }
      }

      if (orderedClusters.size() < 5)
         continue;

      hNClusters->Fill(orderedClusters.size());

      // --- Find MC truth at the first cluster position ---
      // The first cluster is ~10 mm from the vertex. Use the MC point closest
      // to the first cluster as the "local truth" for seeding and comparison.
      XYZPoint firstCluster = orderedClusters.front();
      double bestMCDist = 1e9;
      int bestMCIdx = 0;
      for (size_t j = 0; j < mcProtonPoints.size(); j++) {
         double d = (mcProtonPoints[j].pos - firstCluster).R();
         if (d < bestMCDist) {
            bestMCDist = d;
            bestMCIdx = j;
         }
      }
      XYZVector trueMomLocal = mcProtonPoints[bestMCIdx].mom; // MC momentum near first cluster

      // --- Prepare UKF initial state ---
      double pTrue = trueMomLocal.R();
      double thetaTrue = trueMomLocal.Theta();
      double phiTrue = trueMomLocal.Phi();

      XYZPoint initialPos = firstCluster;
      XYZVector initialMom = trueMomLocal;

      // Initial covariance
      double sigma_pos = 2.0; // mm (wider for digitized data)
      double sigma_mom = 0.15 * pTrue;
      double sigma_ang = 5.0 * M_PI / 180.0; // 5 deg (wider for digi)
      TMatrixD cov(6, 6);
      cov.Zero();
      for (int i = 0; i < 3; i++)
         cov(i, i) = sigma_pos * sigma_pos;
      cov(3, 3) = sigma_mom * sigma_mom;
      cov(4, 4) = sigma_ang * sigma_ang;
      cov(5, 5) = sigma_ang * sigma_ang;

      // Measurement covariance
      double measSigma = 2.0; // mm — fallback for fixed covariance
      TMatrixD defaultMeasCov(3, 3);
      defaultMeasCov.Zero();
      for (int i = 0; i < 3; i++)
         defaultMeasCov(i, i) = measSigma * measSigma;

      // --- Run UKF ---
      ukf.SetInitialState(initialPos, initialMom, cov);
      ukf.SetMeasCov(defaultMeasCov); // Set default (used for first step or if per-cluster disabled)

      bool converged = true;
      try {
         for (size_t i = 1; i < orderedClusters.size(); i++) {
            // Update measurement covariance from cluster if enabled
            if (usePerClusterCov) {
               const TMatrixDSym &clCov = clusters->at(orderedIndices[i]).GetCovMatrix();
               TMatrixD mc(3, 3);
               for (int r = 0; r < 3; r++)
                  for (int c = 0; c < 3; c++)
                     mc(r, c) = clCov(r, c);
               // Floor to avoid singular matrix
               for (int d = 0; d < 3; d++)
                  if (mc(d, d) < 0.01)
                     mc(d, d) = 0.01;
               ukf.SetMeasCov(mc);
            }
            ukf.predictUKF(orderedClusters[i]);
            ukf.correctUKF(orderedClusters[i]);
         }
         ukf.smoothUKF();
      } catch (const std::exception &e) {
         converged = false;
      }

      if (!converged) {
         nFailed++;
         std::cout << "Event " << iEvt << ": FAILED (" << clusters->size() << " clusters)" << std::endl;
         continue;
      }
      nFitted++;

      // --- Extract results ---
      auto &smoothed = ukf.GetSmoothedStates();
      double pReco = smoothed[0][3];
      double thetaReco = smoothed[0][4];

      double momErr = (pReco - pTrue) / pTrue * 100;
      double thetaErr = (thetaReco - thetaTrue) * 180.0 / M_PI;

      // Mean residual between smoothed and measured cluster positions
      double sumResid = 0;
      int nResid = std::min(smoothed.size(), orderedClusters.size()) - 1;
      for (int i = 1; i <= nResid; i++) {
         XYZPoint sp(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
         sumResid += (sp - orderedClusters[i]).R();
      }
      double meanResid = nResid > 0 ? sumResid / nResid : 0;

      hMomErr->Fill(momErr);
      hThetaErr->Fill(thetaErr);
      hResid->Fill(meanResid);

      std::cout << "Event " << iEvt << ": p_true=" << pTrue << " p_reco=" << pReco << " MeV/c"
                << " err=" << momErr << "%" << " theta_err=" << thetaErr << " deg"
                << " resid=" << meanResid << " mm"
                << " clusters=" << clusters->size() << std::endl;
   }

   // --- Summary ---
   std::cout << "\n========================================" << std::endl;
   std::cout << " UKF Digitized Data Validation" << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "  Fitted:  " << nFitted << std::endl;
   std::cout << "  Failed:  " << nFailed << std::endl;
   if (nFitted > 0) {
      std::cout << "  Mom error mean:   " << hMomErr->GetMean() << "%" << std::endl;
      std::cout << "  Mom error RMS:    " << hMomErr->GetRMS() << "%" << std::endl;
      std::cout << "  Theta error mean: " << hThetaErr->GetMean() << " deg" << std::endl;
      std::cout << "  Mean residual:    " << hResid->GetMean() << " mm" << std::endl;
      std::cout << "  Avg clusters:     " << hNClusters->GetMean() << std::endl;
   }
   std::cout << "========================================" << std::endl;

   // --- Save plots ---
   TCanvas *c = new TCanvas("cUKFDigi", "UKF Digitized Validation", 1200, 800);
   c->Divide(2, 2);
   c->cd(1);
   hMomErr->Draw("hist");
   c->cd(2);
   hThetaErr->Draw("hist");
   c->cd(3);
   hResid->Draw("hist");
   c->cd(4);
   hNClusters->Draw("hist");
   c->SaveAs("data/ukf_digi_validation.png");
   std::cout << "Saved: data/ukf_digi_validation.png" << std::endl;
}
