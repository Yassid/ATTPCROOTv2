/// @file display_ukf.C
/// @brief Phase 1: Standalone visualization of UKF-fitted tracks.
///
/// Reads digitized data (output_digi.root) and UKF-fitted output
/// (output_ukf_only.root), then displays clusters vs fitted track
/// with momentum profile and residual plots.
///
/// Prerequisites:
///   1. Run run_digi_attpc.C → output_digi.root
///   2. Run run_ukf_only.C  → output_ukf_only.root
///
/// Run: root -l display_ukf.C
///      root -l 'display_ukf.C(3)' // display event 3

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

const double ZPadPlane = 1000.0;

// Global state for event navigation
TFile *gFileDigi = nullptr;
TFile *gFileFit = nullptr;
TTree *gTreeDigi = nullptr;
TTree *gTreeFit = nullptr;
TClonesArray *gPatEvtArr = nullptr;
TClonesArray *gTrackingEvtArr = nullptr;
int gCurrentEvent = 0;
int gNEvents = 0;

// Canvases
TCanvas *gCanvas3D = nullptr;
TCanvas *gCanvasDiag = nullptr;

void DrawEvent(int iEvt);

// ---------------------------------------------------------------------------
void display_ukf(int startEvent = 1)
{
   TString inDir = "./data/";

   gFileDigi = TFile::Open(inDir + "output_digi.root");
   gFileFit = TFile::Open(inDir + "output_ukf_only.root");

   if (!gFileDigi || !gFileFit) {
      std::cerr << "Missing input files. Run run_digi_attpc.C and run_ukf_only.C first." << std::endl;
      return;
   }

   gTreeDigi = (TTree *)gFileDigi->Get("cbmsim");
   gTreeFit = (TTree *)gFileFit->Get("cbmsim");

   gTreeDigi->SetBranchAddress("AtPatternEvent", &gPatEvtArr);
   gTreeFit->SetBranchAddress("AtTrackingEvent", &gTrackingEvtArr);

   gNEvents = gTreeDigi->GetEntries();
   gCurrentEvent = startEvent;

   // Create canvases
   gCanvas3D = new TCanvas("c3D", "UKF Track Display", 1200, 600);
   gCanvas3D->Divide(3, 1);

   gCanvasDiag = new TCanvas("cDiag", "UKF Diagnostics", 1000, 600);
   gCanvasDiag->Divide(2, 2);

   DrawEvent(gCurrentEvent);

   std::cout << "\n=== UKF Track Display ===" << std::endl;
   std::cout << "Events: " << gNEvents << std::endl;
   std::cout << "Commands:" << std::endl;
   std::cout << "  DrawEvent(n)  — go to event n" << std::endl;
   std::cout << "  DrawEvent(-1) — next event with a track" << std::endl;
   std::cout << "  DrawEvent(-2) — previous event with a track" << std::endl;
}

// ---------------------------------------------------------------------------
void DrawEvent(int iEvt)
{
   // Navigation
   if (iEvt == -1) { // next with track
      for (int i = gCurrentEvent + 1; i < gNEvents; i++) {
         gTreeDigi->GetEntry(i);
         if (gPatEvtArr->GetEntries() > 0) {
            AtPatternEvent *pe = (AtPatternEvent *)gPatEvtArr->At(0);
            if (!pe->GetTrackCand().empty() && pe->GetTrackCand()[0].GetHitClusterArray()->size() > 10) {
               iEvt = i;
               break;
            }
         }
      }
      if (iEvt == -1) {
         std::cout << "No more events with tracks." << std::endl;
         return;
      }
   } else if (iEvt == -2) { // prev with track
      for (int i = gCurrentEvent - 1; i >= 0; i--) {
         gTreeDigi->GetEntry(i);
         if (gPatEvtArr->GetEntries() > 0) {
            AtPatternEvent *pe = (AtPatternEvent *)gPatEvtArr->At(0);
            if (!pe->GetTrackCand().empty() && pe->GetTrackCand()[0].GetHitClusterArray()->size() > 10) {
               iEvt = i;
               break;
            }
         }
      }
      if (iEvt == -2) {
         std::cout << "No previous events with tracks." << std::endl;
         return;
      }
   }

   gCurrentEvent = iEvt;
   gTreeDigi->GetEntry(iEvt);
   gTreeFit->GetEntry(iEvt);

   // Get pattern event
   if (gPatEvtArr->GetEntries() == 0) {
      std::cout << "Event " << iEvt << ": no data" << std::endl;
      return;
   }
   AtPatternEvent *patEvt = (AtPatternEvent *)gPatEvtArr->At(0);
   std::vector<AtTrack> &tracks = patEvt->GetTrackCand();

   // Get fitted tracks
   AtTrackingEvent *trackingEvt = nullptr;
   const std::vector<std::unique_ptr<AtFittedTrack>> *fittedTracks = nullptr;
   if (gTrackingEvtArr && gTrackingEvtArr->GetEntries() > 0) {
      trackingEvt = (AtTrackingEvent *)gTrackingEvtArr->At(0);
      fittedTracks = &trackingEvt->GetFittedTracks();
   }

   if (tracks.empty()) {
      std::cout << "Event " << iEvt << ": no tracks" << std::endl;
      return;
   }

   // Use first (largest) track
   AtTrack *track = &tracks[0];
   std::vector<AtHitCluster> *clusters = track->GetHitClusterArray();

   // Convert clusters to lab frame
   std::vector<double> cx, cy, cz;
   for (auto &cl : *clusters) {
      XYZPoint pos = cl.GetPosition();
      cx.push_back(pos.X());
      cy.push_back(pos.Y());
      cz.push_back(ZPadPlane - pos.Z()); // digi → lab
   }

   // Get smoothed positions from fitted track
   std::vector<double> sx, sy, sz;
   std::vector<double> momProfile;
   std::vector<double> residuals;
   double pReco = -1, thetaReco = -1, phiReco = -1, keReco = -1;
   bool hasFit = false;

   if (fittedTracks && !fittedTracks->empty()) {
      AtFittedTrack *ft = fittedTracks->at(0).get();
      auto &smoothed = ft->GetSmoothedPositions();
      auto kin = ft->GetKinematics();
      pReco = std::sqrt(2 * 938.272 * kin.kineticEnergy + kin.kineticEnergy * kin.kineticEnergy);
      thetaReco = kin.theta * 180.0 / M_PI;
      phiReco = kin.phi * 180.0 / M_PI;
      keReco = kin.kineticEnergy;

      // Smoothed positions are in lab frame (after Z conversion in AtFitterUKF)
      XYZVector vtx = ft->GetVertex();
      sx.push_back(vtx.X());
      sy.push_back(vtx.Y());
      sz.push_back(vtx.Z());

      for (auto &sp : smoothed) {
         sx.push_back(sp.X());
         sy.push_back(sp.Y());
         sz.push_back(sp.Z());
      }

      // Residuals (smoothed vs measured clusters)
      int nPts = std::min(sx.size(), cx.size());
      for (int i = 0; i < nPts; i++) {
         double dx = sx[i] - cx[i];
         double dy = sy[i] - cy[i];
         double dz = sz[i] - cz[i];
         residuals.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
      }

      hasFit = true;
   }

   // ===== 3D Canvas (2D projections) =====
   gCanvas3D->cd(1);
   TGraph *gXY_cl = new TGraph(cx.size(), cx.data(), cy.data());
   gXY_cl->SetTitle(Form("Event %d  XY;X [mm];Y [mm]", iEvt));
   gXY_cl->SetMarkerStyle(20);
   gXY_cl->SetMarkerSize(0.6);
   gXY_cl->SetMarkerColor(kBlue);
   gXY_cl->Draw("AP");
   if (hasFit) {
      TGraph *gXY_fit = new TGraph(sx.size(), sx.data(), sy.data());
      gXY_fit->SetMarkerStyle(21);
      gXY_fit->SetMarkerSize(0.5);
      gXY_fit->SetMarkerColor(kRed);
      gXY_fit->Draw("P SAME");
      TLegend *leg1 = new TLegend(0.15, 0.75, 0.45, 0.88);
      leg1->AddEntry(gXY_cl, "Clusters", "p");
      leg1->AddEntry(gXY_fit, "Smoothed", "p");
      leg1->Draw();
   }

   gCanvas3D->cd(2);
   TGraph *gXZ_cl = new TGraph(cx.size(), cz.data(), cx.data());
   gXZ_cl->SetTitle("XZ;Z_{lab} [mm];X [mm]");
   gXZ_cl->SetMarkerStyle(20);
   gXZ_cl->SetMarkerSize(0.6);
   gXZ_cl->SetMarkerColor(kBlue);
   gXZ_cl->Draw("AP");
   if (hasFit) {
      TGraph *gXZ_fit = new TGraph(sx.size(), sz.data(), sx.data());
      gXZ_fit->SetMarkerStyle(21);
      gXZ_fit->SetMarkerSize(0.5);
      gXZ_fit->SetMarkerColor(kRed);
      gXZ_fit->Draw("P SAME");
   }

   gCanvas3D->cd(3);
   TGraph *gYZ_cl = new TGraph(cy.size(), cz.data(), cy.data());
   gYZ_cl->SetTitle("YZ;Z_{lab} [mm];Y [mm]");
   gYZ_cl->SetMarkerStyle(20);
   gYZ_cl->SetMarkerSize(0.6);
   gYZ_cl->SetMarkerColor(kBlue);
   gYZ_cl->Draw("AP");
   if (hasFit) {
      TGraph *gYZ_fit = new TGraph(sy.size(), sz.data(), sy.data());
      gYZ_fit->SetMarkerStyle(21);
      gYZ_fit->SetMarkerSize(0.5);
      gYZ_fit->SetMarkerColor(kRed);
      gYZ_fit->Draw("P SAME");
   }

   gCanvas3D->Update();

   // ===== Diagnostics Canvas =====
   if (hasFit) {
      // Residuals vs cluster index
      gCanvasDiag->cd(1);
      std::vector<double> indices(residuals.size());
      for (int i = 0; i < (int)indices.size(); i++)
         indices[i] = i;
      TGraph *gResid = new TGraph(residuals.size(), indices.data(), residuals.data());
      gResid->SetTitle("Residuals;Cluster index;|smoothed - measured| [mm]");
      gResid->SetMarkerStyle(20);
      gResid->SetMarkerSize(0.5);
      gResid->SetMarkerColor(kMagenta);
      gResid->Draw("AP");

      // Distance from vertex along track
      gCanvasDiag->cd(2);
      std::vector<double> arcLen(cx.size());
      arcLen[0] = 0;
      for (int i = 1; i < (int)cx.size(); i++) {
         double dx = cx[i] - cx[i - 1];
         double dy = cy[i] - cy[i - 1];
         double dz = cz[i] - cz[i - 1];
         arcLen[i] = arcLen[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
      }
      TGraph *gArc = new TGraph(cx.size(), arcLen.data(), cz.data());
      gArc->SetTitle("Track arc length;Arc length [mm];Z_{lab} [mm]");
      gArc->SetMarkerStyle(20);
      gArc->SetMarkerSize(0.5);
      gArc->Draw("AP");

      // Cluster charge profile
      gCanvasDiag->cd(3);
      std::vector<double> charges(clusters->size());
      for (int i = 0; i < (int)clusters->size(); i++)
         charges[i] = clusters->at(i).GetCharge();
      TGraph *gQ = new TGraph(charges.size(), indices.data(), charges.data());
      gQ->SetTitle("Charge profile;Cluster index;Charge [arb]");
      gQ->SetMarkerStyle(20);
      gQ->SetMarkerSize(0.5);
      gQ->SetMarkerColor(kGreen + 2);
      gQ->Draw("AP");

      // Fit info text
      gCanvasDiag->cd(4);
      TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
      pt->SetTextAlign(12);
      pt->SetTextSize(0.06);
      pt->SetFillColor(0);
      pt->AddText(Form("Event %d", iEvt));
      pt->AddText(Form("Clusters: %d", (int)clusters->size()));
      pt->AddText(Form("p = %.1f MeV/c", pReco));
      pt->AddText(Form("KE = %.2f MeV", keReco));
      pt->AddText(Form("#theta = %.1f deg", thetaReco));
      pt->AddText(Form("#phi = %.1f deg", phiReco));
      double meanResid = 0;
      for (double r : residuals)
         meanResid += r;
      if (!residuals.empty())
         meanResid /= residuals.size();
      pt->AddText(Form("Mean residual: %.2f mm", meanResid));
      pt->Draw();
   } else {
      gCanvasDiag->cd(4);
      TPaveText *pt = new TPaveText(0.1, 0.3, 0.9, 0.7, "NDC");
      pt->SetTextSize(0.08);
      pt->SetFillColor(0);
      pt->AddText(Form("Event %d", iEvt));
      pt->AddText(Form("Clusters: %d", (int)clusters->size()));
      pt->AddText("No fit result");
      pt->Draw();
   }

   gCanvasDiag->Update();

   std::cout << "Event " << iEvt << ": " << clusters->size() << " clusters";
   if (hasFit)
      std::cout << ", p=" << pReco << " MeV/c, KE=" << keReco << " MeV, theta=" << thetaReco << " deg";
   std::cout << std::endl;
}
