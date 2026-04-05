/// @file display_ukf_gnn.C
/// @brief Interactive display of GNN-clustered UKF-fitted tracks.
///
/// Same interface as display_ukf.C but reads GNN clusters from CSV
/// and GNN fit results from output_ukf_gnn.root.
///
/// Prerequisites:
///   1. Run gnn_cluster_export.py → gnn_clusters.csv
///   2. Run run_ukf_gnn.C → output_ukf_gnn.root
///
/// Run: root -l display_ukf_gnn.C
///      root -l 'display_ukf_gnn.C(3)' // display event 3

#include <fstream>
#include <sstream>
#include <map>

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

const double ZPadPlane = 1000.0;

// ─── GNN cluster data ──────────────────────────────────────────────

struct GNNCluster {
   int cluster_id, track_pred, n_hits;
   double x, y, z, charge;
};

std::map<int, std::vector<GNNCluster>> gGNNClusters;

void LoadGNNClusters()
{
   std::ifstream f("data/gnn_training/gnn_clusters.csv");
   if (!f.is_open()) { std::cerr << "Cannot open gnn_clusters.csv" << std::endl; return; }
   std::string line;
   std::getline(f, line);
   while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string tok;
      int eid;
      std::getline(ss, tok, ','); eid = std::stoi(tok);
      GNNCluster cl;
      std::getline(ss, tok, ','); cl.cluster_id = std::stoi(tok);
      std::getline(ss, tok, ','); cl.track_pred = std::stoi(tok);
      std::getline(ss, tok, ','); cl.n_hits = std::stoi(tok);
      std::getline(ss, tok, ','); cl.x = std::stod(tok);
      std::getline(ss, tok, ','); cl.y = std::stod(tok);
      std::getline(ss, tok, ','); cl.z = std::stod(tok);
      std::getline(ss, tok, ','); cl.charge = std::stod(tok);
      // skip covariance columns
      gGNNClusters[eid].push_back(cl);
   }
   std::cout << "Loaded GNN clusters for " << gGNNClusters.size() << " events" << std::endl;
}

// ─── GNN fit results indexed by event_id ───────────────────────────

struct GNNFitResult {
   double KE, theta, phi, vtxX, vtxY, vtxZ, length;
   int nClusters, status;
};

std::map<int, GNNFitResult> gGNNFits;

void LoadGNNFits()
{
   TFile *f = TFile::Open("data/output_ukf_gnn.root");
   if (!f) { std::cerr << "Cannot open output_ukf_gnn.root" << std::endl; return; }
   TTree *t = (TTree *)f->Get("convergence");

   double KE, theta, phi, vtxX, vtxY, vtxZ, length;
   int nCl, status, evid, trackPred;
   t->SetBranchAddress("event_id", &evid);
   t->SetBranchAddress("KE", &KE);
   t->SetBranchAddress("theta", &theta);
   t->SetBranchAddress("phi", &phi);
   t->SetBranchAddress("vtxX", &vtxX);
   t->SetBranchAddress("vtxY", &vtxY);
   t->SetBranchAddress("vtxZ", &vtxZ);
   t->SetBranchAddress("length", &length);
   t->SetBranchAddress("nClusters", &nCl);
   t->SetBranchAddress("status", &status);

   for (int i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (status == 0 && KE > 0) {
         gGNNFits[evid] = {KE, theta, phi, vtxX, vtxY, vtxZ, length, nCl, status};
      }
   }
   f->Close();
   std::cout << "Loaded GNN fit results for " << gGNNFits.size() << " events" << std::endl;
}

// ─── Global state ──────────────────────────────────────────────────

int gCurrentEvent = 0;
int gNEvents = 0;
TCanvas *gCanvas3D = nullptr;
TCanvas *gCanvasDiag = nullptr;

// List of events that have GNN fits (for navigation)
std::vector<int> gFittedEvents;
int gFittedIdx = 0;

void DrawEvent(int iEvt);

// ─── Init ──────────────────────────────────────────────────────────

void display_ukf_gnn(int startEvent = -1)
{
   LoadGNNClusters();
   LoadGNNFits();

   // Build sorted list of fitted events
   for (auto &[eid, fit] : gGNNFits)
      gFittedEvents.push_back(eid);
   std::sort(gFittedEvents.begin(), gFittedEvents.end());
   gNEvents = gFittedEvents.size();

   gCanvas3D = new TCanvas("c3D", "GNN UKF Track Display", 1200, 600);
   gCanvas3D->Divide(3, 1);

   gCanvasDiag = new TCanvas("cDiag", "GNN UKF Diagnostics", 1000, 600);
   gCanvasDiag->Divide(2, 2);

   if (startEvent < 0 && !gFittedEvents.empty())
      startEvent = gFittedEvents[0];

   DrawEvent(startEvent);

   std::cout << "\n=== GNN UKF Track Display ===" << std::endl;
   std::cout << "Fitted events: " << gNEvents << std::endl;
   std::cout << "Commands:" << std::endl;
   std::cout << "  DrawEvent(n)  — go to event n" << std::endl;
   std::cout << "  DrawEvent(-1) — next fitted event" << std::endl;
   std::cout << "  DrawEvent(-2) — previous fitted event" << std::endl;
}

// ─── Draw ──────────────────────────────────────────────────────────

void DrawEvent(int iEvt)
{
   // Navigation: cycle through fitted events only
   if (iEvt == -1) {
      gFittedIdx++;
      if (gFittedIdx >= (int)gFittedEvents.size()) {
         std::cout << "No more fitted events." << std::endl;
         gFittedIdx = gFittedEvents.size() - 1;
         return;
      }
      iEvt = gFittedEvents[gFittedIdx];
   } else if (iEvt == -2) {
      gFittedIdx--;
      if (gFittedIdx < 0) {
         std::cout << "No previous fitted events." << std::endl;
         gFittedIdx = 0;
         return;
      }
      iEvt = gFittedEvents[gFittedIdx];
   } else {
      // Find this event in the fitted list
      auto it = std::lower_bound(gFittedEvents.begin(), gFittedEvents.end(), iEvt);
      if (it == gFittedEvents.end() || *it != iEvt) {
         // Find nearest fitted event
         if (it != gFittedEvents.end()) {
            iEvt = *it;
            gFittedIdx = it - gFittedEvents.begin();
         } else if (!gFittedEvents.empty()) {
            iEvt = gFittedEvents.back();
            gFittedIdx = gFittedEvents.size() - 1;
         } else {
            std::cout << "No fitted events." << std::endl;
            return;
         }
         std::cout << "Jumped to nearest fitted event: " << iEvt << std::endl;
      } else {
         gFittedIdx = it - gFittedEvents.begin();
      }
   }

   gCurrentEvent = iEvt;

   // Get GNN clusters for this event — find the track with most clusters
   auto itCl = gGNNClusters.find(iEvt);
   if (itCl == gGNNClusters.end()) {
      std::cout << "Event " << iEvt << ": no GNN clusters" << std::endl;
      return;
   }

   // Find the largest track
   std::map<int, int> trackSizes;
   for (auto &cl : itCl->second) trackSizes[cl.track_pred]++;
   int bestTrack = -1, maxSize = 0;
   for (auto &[tp, sz] : trackSizes)
      if (sz > maxSize) { maxSize = sz; bestTrack = tp; }

   // Collect clusters for the best track, convert to lab frame
   std::vector<double> cx, cy, cz, charges;
   for (auto &cl : itCl->second) {
      if (cl.track_pred == bestTrack) {
         cx.push_back(cl.x);
         cy.push_back(cl.y);
         cz.push_back(ZPadPlane - cl.z);  // digi → lab
         charges.push_back(cl.charge);
      }
   }

   // Get fit result
   auto itFit = gGNNFits.find(iEvt);
   bool hasFit = (itFit != gGNNFits.end());
   double pReco = -1, thetaReco = -1, phiReco = -1, keReco = -1, trackLen = -1;

   if (hasFit) {
      auto &fit = itFit->second;
      keReco = fit.KE;
      thetaReco = fit.theta * 180.0 / M_PI;
      phiReco = fit.phi * 180.0 / M_PI;
      pReco = std::sqrt(2 * 938.272 * keReco + keReco * keReco);
      trackLen = fit.length;
   }

   // ═══════ Track Canvas (XY, XZ, YZ) ═══════
   gCanvas3D->cd(1);
   gPad->Clear();
   TGraph *gXY = new TGraph(cx.size(), cx.data(), cy.data());
   gXY->SetTitle(Form("Event %d — GNN XY;X [mm];Y [mm]", iEvt));
   gXY->SetMarkerStyle(20); gXY->SetMarkerSize(0.6); gXY->SetMarkerColor(kGreen + 2);
   gXY->Draw("AP");

   gCanvas3D->cd(2);
   gPad->Clear();
   TGraph *gXZ = new TGraph(cx.size(), cz.data(), cx.data());
   gXZ->SetTitle("GNN ZX;Z_{lab} [mm];X [mm]");
   gXZ->SetMarkerStyle(20); gXZ->SetMarkerSize(0.6); gXZ->SetMarkerColor(kGreen + 2);
   gXZ->Draw("AP");

   gCanvas3D->cd(3);
   gPad->Clear();
   TGraph *gYZ = new TGraph(cy.size(), cz.data(), cy.data());
   gYZ->SetTitle("GNN ZY;Z_{lab} [mm];Y [mm]");
   gYZ->SetMarkerStyle(20); gYZ->SetMarkerSize(0.6); gYZ->SetMarkerColor(kGreen + 2);
   gYZ->Draw("AP");

   gCanvas3D->Update();

   // ═══════ Diagnostics Canvas ═══════
   // Charge profile
   gCanvasDiag->cd(1);
   gPad->Clear();
   std::vector<double> indices(charges.size());
   for (int i = 0; i < (int)indices.size(); i++) indices[i] = i;
   TGraph *gQ = new TGraph(charges.size(), indices.data(), charges.data());
   gQ->SetTitle("Charge profile;Cluster index;Charge [arb]");
   gQ->SetMarkerStyle(20); gQ->SetMarkerSize(0.5); gQ->SetMarkerColor(kGreen + 2);
   gQ->Draw("AP");

   // Z profile along arc
   gCanvasDiag->cd(2);
   gPad->Clear();
   std::vector<double> arcLen(cx.size());
   arcLen[0] = 0;
   for (int i = 1; i < (int)cx.size(); i++) {
      double dx = cx[i] - cx[i - 1];
      double dy = cy[i] - cy[i - 1];
      double dz = cz[i] - cz[i - 1];
      arcLen[i] = arcLen[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
   }
   TGraph *gArc = new TGraph(cx.size(), arcLen.data(), cz.data());
   gArc->SetTitle("Track arc;Arc length [mm];Z_{lab} [mm]");
   gArc->SetMarkerStyle(20); gArc->SetMarkerSize(0.5);
   gArc->Draw("AP");

   // Inter-cluster spacing
   gCanvasDiag->cd(3);
   gPad->Clear();
   if (cx.size() > 1) {
      std::vector<double> spacing(cx.size() - 1);
      std::vector<double> spIdx(cx.size() - 1);
      for (int i = 0; i < (int)spacing.size(); i++) {
         double dx = cx[i + 1] - cx[i];
         double dy = cy[i + 1] - cy[i];
         double dz = cz[i + 1] - cz[i];
         spacing[i] = std::sqrt(dx * dx + dy * dy + dz * dz);
         spIdx[i] = i;
      }
      TGraph *gSp = new TGraph(spacing.size(), spIdx.data(), spacing.data());
      gSp->SetTitle("Inter-cluster spacing;Cluster index;Spacing [mm]");
      gSp->SetMarkerStyle(20); gSp->SetMarkerSize(0.5); gSp->SetMarkerColor(kMagenta);
      gSp->Draw("AP");
   }

   // Fit info
   gCanvasDiag->cd(4);
   gPad->Clear();
   TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
   pt->SetTextAlign(12); pt->SetTextSize(0.06); pt->SetFillColor(0);
   pt->AddText(Form("#bf{Event %d} (GNN Arc-Walk)", iEvt));
   pt->AddText(Form("Clusters: %d  (%d hits/cluster)", (int)cx.size(), cx.size() > 0 ? 8 : 0));
   pt->AddText(Form("Arc length: %.0f mm", arcLen.empty() ? 0 : arcLen.back()));
   pt->AddText("");
   if (hasFit) {
      pt->AddText(Form("p = %.1f MeV/c", pReco));
      pt->AddText(Form("KE = %.2f MeV", keReco));
      pt->AddText(Form("#theta = %.1f#circ", thetaReco));
      pt->AddText(Form("#phi = %.1f#circ", phiReco));
      pt->AddText(Form("Track length = %.0f mm", trackLen));
   } else {
      pt->AddText("No fit result");
   }
   pt->Draw();

   gCanvasDiag->Update();

   std::cout << "Event " << iEvt << " [" << gFittedIdx + 1 << "/" << gNEvents << "]: "
             << cx.size() << " clusters";
   if (hasFit)
      std::cout << ", p=" << Form("%.1f", pReco) << " MeV/c, KE=" << Form("%.2f", keReco)
                << " MeV, theta=" << Form("%.1f", thetaReco) << " deg";
   std::cout << std::endl;
}
