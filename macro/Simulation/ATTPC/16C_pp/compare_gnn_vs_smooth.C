/// @file compare_gnn_vs_smooth.C
/// @brief Compare UKF fit results: GNN arc-walk clusters vs ClusterizeSmooth3D.
///
/// Reads output_ukf_gnn.root (flat TTree from run_ukf_gnn.C) and
/// output_ukf_only.root (AtTrackingEvent from run_ukf_only.C),
/// overlays KE, theta, vertex distributions.
///
/// Run: root -b -q compare_gnn_vs_smooth.C

void compare_gnn_vs_smooth()
{
   // ── 1. Load GNN results (flat tree) ──
   TFile *fGNN = TFile::Open("data/output_ukf_gnn.root");
   if (!fGNN) { std::cerr << "Cannot open output_ukf_gnn.root" << std::endl; return; }
   TTree *tGNN = (TTree *)fGNN->Get("convergence");

   double gnn_KE, gnn_theta, gnn_phi, gnn_vtxX, gnn_vtxY, gnn_vtxZ, gnn_length;
   int gnn_nCl, gnn_status, gnn_evid;
   tGNN->SetBranchAddress("KE", &gnn_KE);
   tGNN->SetBranchAddress("theta", &gnn_theta);
   tGNN->SetBranchAddress("phi", &gnn_phi);
   tGNN->SetBranchAddress("vtxX", &gnn_vtxX);
   tGNN->SetBranchAddress("vtxY", &gnn_vtxY);
   tGNN->SetBranchAddress("vtxZ", &gnn_vtxZ);
   tGNN->SetBranchAddress("length", &gnn_length);
   tGNN->SetBranchAddress("nClusters", &gnn_nCl);
   tGNN->SetBranchAddress("status", &gnn_status);
   tGNN->SetBranchAddress("event_id", &gnn_evid);

   // ── 2. Load ClusterizeSmooth3D results ──
   TFile *fSmooth = TFile::Open("data/output_ukf_only.root");
   if (!fSmooth) { std::cerr << "Cannot open output_ukf_only.root" << std::endl; return; }
   TTree *tSmooth = (TTree *)fSmooth->Get("cbmsim");

   TClonesArray *trackEvtArr = nullptr;
   tSmooth->SetBranchAddress("AtTrackingEvent", &trackEvtArr);

   // ── 3. Collect data ──
   std::vector<double> sKE, sTheta, sPhi, sVtxR;
   std::vector<double> gKE, gTheta, gPhi, gVtxR;

   // GNN data
   for (int i = 0; i < tGNN->GetEntries(); i++) {
      tGNN->GetEntry(i);
      if (gnn_status != 0 || gnn_KE < 0) continue;
      gKE.push_back(gnn_KE);
      gTheta.push_back(gnn_theta * 180.0 / M_PI); // rad → deg
      gPhi.push_back(gnn_phi * 180.0 / M_PI);
      gVtxR.push_back(std::sqrt(gnn_vtxX * gnn_vtxX + gnn_vtxY * gnn_vtxY));
   }

   // Smooth data
   for (int i = 0; i < tSmooth->GetEntries(); i++) {
      tSmooth->GetEntry(i);
      if (!trackEvtArr || trackEvtArr->GetEntries() == 0) continue;
      auto *trackEvt = (AtTrackingEvent *)trackEvtArr->At(0);
      auto &fittedTracks = trackEvt->GetFittedTracks();
      for (auto &ft : fittedTracks) {
         auto &kin = ft->GetKinematics();
         if (kin.kineticEnergy < 0) continue;
         sKE.push_back(kin.kineticEnergy);
         sTheta.push_back(kin.theta * 180.0 / M_PI);
         sPhi.push_back(kin.phi * 180.0 / M_PI);
         auto vtx = ft->GetVertex();
         sVtxR.push_back(std::sqrt(vtx.X() * vtx.X() + vtx.Y() * vtx.Y()));
      }
   }

   std::cout << "Smooth3D tracks: " << sKE.size() << std::endl;
   std::cout << "GNN tracks:      " << gKE.size() << std::endl;

   // ── 4. Plot ──
   auto *c1 = new TCanvas("c1", "GNN vs Smooth3D Comparison", 1600, 1000);
   c1->Divide(3, 2);

   // KE distributions
   c1->cd(1);
   auto *hSKE = new TH1D("hSKE", "Kinetic Energy;KE [MeV];Counts", 100, 0, 15);
   auto *hGKE = new TH1D("hGKE", "Kinetic Energy;KE [MeV];Counts", 100, 0, 15);
   for (auto v : sKE) hSKE->Fill(v);
   for (auto v : gKE) hGKE->Fill(v);
   hSKE->SetLineColor(kBlue);
   hGKE->SetLineColor(kRed);
   hSKE->SetLineWidth(2);
   hGKE->SetLineWidth(2);
   hSKE->Draw();
   hGKE->Draw("SAME");
   auto *leg1 = new TLegend(0.55, 0.7, 0.9, 0.9);
   leg1->AddEntry(hSKE, Form("Smooth3D (N=%lu)", sKE.size()), "l");
   leg1->AddEntry(hGKE, Form("GNN (N=%lu)", gKE.size()), "l");
   leg1->Draw();

   // Theta distributions
   c1->cd(2);
   auto *hST = new TH1D("hST", "Scattering Angle;#theta [deg];Counts", 100, 0, 180);
   auto *hGT = new TH1D("hGT", "Scattering Angle;#theta [deg];Counts", 100, 0, 180);
   for (auto v : sTheta) hST->Fill(v);
   for (auto v : gTheta) hGT->Fill(v);
   hST->SetLineColor(kBlue);
   hGT->SetLineColor(kRed);
   hST->SetLineWidth(2);
   hGT->SetLineWidth(2);
   hST->Draw();
   hGT->Draw("SAME");
   auto *leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
   leg2->AddEntry(hST, "Smooth3D", "l");
   leg2->AddEntry(hGT, "GNN", "l");
   leg2->Draw();

   // Vertex R distributions
   c1->cd(3);
   auto *hSV = new TH1D("hSV", "Vertex R_{xy};R [mm];Counts", 100, 0, 50);
   auto *hGV = new TH1D("hGV", "Vertex R_{xy};R [mm];Counts", 100, 0, 50);
   for (auto v : sVtxR) hSV->Fill(v);
   for (auto v : gVtxR) hGV->Fill(v);
   hSV->SetLineColor(kBlue);
   hGV->SetLineColor(kRed);
   hSV->SetLineWidth(2);
   hGV->SetLineWidth(2);
   hSV->Draw();
   hGV->Draw("SAME");
   auto *leg3 = new TLegend(0.55, 0.7, 0.9, 0.9);
   leg3->AddEntry(hSV, "Smooth3D", "l");
   leg3->AddEntry(hGV, "GNN", "l");
   leg3->Draw();

   // KE vs Theta (kinematic curve)
   c1->cd(4);
   auto *hSKT = new TH2D("hSKT", "Smooth3D: KE vs #theta;#theta [deg];KE [MeV]", 90, 0, 180, 75, 0, 15);
   for (size_t i = 0; i < sKE.size(); i++) hSKT->Fill(sTheta[i], sKE[i]);
   hSKT->Draw("COLZ");

   c1->cd(5);
   auto *hGKT = new TH2D("hGKT", "GNN: KE vs #theta;#theta [deg];KE [MeV]", 90, 0, 180, 75, 0, 15);
   for (size_t i = 0; i < gKE.size(); i++) hGKT->Fill(gTheta[i], gKE[i]);
   hGKT->Draw("COLZ");

   // Track length comparison
   c1->cd(6);
   auto *hGL = new TH1D("hGL", "GNN Track Length;Length [mm];Counts", 100, 0, 1000);
   for (int i = 0; i < tGNN->GetEntries(); i++) {
      tGNN->GetEntry(i);
      if (gnn_status == 0 && gnn_length > 0) hGL->Fill(gnn_length);
   }
   hGL->SetLineColor(kRed);
   hGL->SetLineWidth(2);
   hGL->Draw();

   c1->SaveAs("data/gnn_vs_smooth_comparison.png");
   std::cout << "Saved: data/gnn_vs_smooth_comparison.png" << std::endl;

   // ── 5. Summary statistics ──
   std::cout << "\n=== Summary ===" << std::endl;
   std::cout << "              Smooth3D       GNN" << std::endl;
   std::cout << Form("Tracks:       %6lu         %6lu", sKE.size(), gKE.size()) << std::endl;
   std::cout << Form("KE mean:      %6.2f MeV    %6.2f MeV", hSKE->GetMean(), hGKE->GetMean()) << std::endl;
   std::cout << Form("KE RMS:       %6.2f MeV    %6.2f MeV", hSKE->GetRMS(), hGKE->GetRMS()) << std::endl;
   std::cout << Form("Theta mean:   %6.1f deg    %6.1f deg", hST->GetMean(), hGT->GetMean()) << std::endl;
   std::cout << Form("VtxR mean:    %6.1f mm     %6.1f mm", hSV->GetMean(), hGV->GetMean()) << std::endl;

   fGNN->Close();
   fSmooth->Close();
}
