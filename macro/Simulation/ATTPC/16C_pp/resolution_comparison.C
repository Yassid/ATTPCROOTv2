/// @file resolution_comparison.C
/// @brief Compare KE and theta resolution for Smooth3D vs GNN+Smooth3D.
///
/// Matches fitted tracks to MC truth proton by closest theta.
/// Uses output_digi.root for MC truth (hit-level trackID + energy).
///
/// Run: root -b -q resolution_comparison.C

void resolution_comparison()
{
   TFile *fDigi = TFile::Open("data/output_digi.root");
   TFile *fSmooth = TFile::Open("data/smooth3d/output_ukf_only.root");
   TFile *fGNN = TFile::Open("data/output_ukf_gnn_smooth.root");

   if (!fDigi || !fSmooth || !fGNN) {
      std::cerr << "Missing files" << std::endl;
      return;
   }

   TTree *tDigi = (TTree *)fDigi->Get("cbmsim");
   TTree *tSmooth = (TTree *)fSmooth->Get("cbmsim");
   TTree *tGNN = (TTree *)fGNN->Get("cbmsim");

   TClonesArray *evtArr = nullptr, *sArr = nullptr, *gArr = nullptr;
   tDigi->SetBranchAddress("AtEventH", &evtArr);
   tSmooth->SetBranchAddress("AtTrackingEvent", &sArr);
   tGNN->SetBranchAddress("AtTrackingEvent", &gArr);

   int nEntries = std::min({(int)tDigi->GetEntries(),
                            (int)tSmooth->GetEntries(),
                            (int)tGNN->GetEntries()});

   // Histograms
   auto *hS_dKE = new TH1D("hS_dKE", "Smooth3D: #DeltaKE/KE;(KE_{fit} - KE_{truth}) / KE_{truth} [%];Counts", 100, -50, 50);
   auto *hG_dKE = new TH1D("hG_dKE", "GNN: #DeltaKE/KE;(KE_{fit} - KE_{truth}) / KE_{truth} [%];Counts", 100, -50, 50);
   auto *hS_dTh = new TH1D("hS_dTh", "Smooth3D: #Delta#theta;#theta_{fit} - #theta_{truth} [deg];Counts", 100, -10, 10);
   auto *hG_dTh = new TH1D("hG_dTh", "GNN: #Delta#theta;#theta_{fit} - #theta_{truth} [deg];Counts", 100, -10, 10);

   // 2D: truth KE vs resolution
   auto *hS_KEvRes = new TH2D("hS_KEvRes", "Smooth3D;KE_{truth} [MeV];#DeltaKE/KE [%]", 50, 0, 15, 100, -50, 50);
   auto *hG_KEvRes = new TH2D("hG_KEvRes", "GNN;KE_{truth} [MeV];#DeltaKE/KE [%]", 50, 0, 15, 100, -50, 50);

   // Per-bin resolution
   auto *hS_KE_vs_truth = new TH2D("hS_KvT", "Smooth3D;KE_{truth} [MeV];KE_{fit} [MeV]", 75, 0, 15, 75, 0, 15);
   auto *hG_KE_vs_truth = new TH2D("hG_KvT", "GNN;KE_{truth} [MeV];KE_{fit} [MeV]", 75, 0, 15, 75, 0, 15);

   int nS = 0, nG = 0, nTruth = 0;

   for (int iEv = 0; iEv < nEntries; iEv++) {
      tDigi->GetEntry(iEv);

      if (!evtArr || evtArr->GetEntries() == 0) continue;
      auto *evt = (AtEvent *)evtArr->At(0);

      // Find proton MC truth: trackID == 1
      double truthKE = -1;
      double truthTheta = -1;
      int nProtonHits = 0;

      for (int h = 0; h < evt->GetNumHits(); h++) {
         auto &hit = evt->GetHit(h);
         auto &mc = hit.GetMCSimPointArray();
         if (!mc.empty() && mc[0].trackID == 1) {
            if (truthKE < 0 && mc[0].energy > 0)
               truthKE = mc[0].energy; // MeV
            nProtonHits++;
         }
      }

      if (truthKE <= 0 || nProtonHits < 10) continue;

      // Compute truth theta from proton hits (first few hits, Z-sorted)
      std::vector<std::pair<double, ROOT::Math::XYZPoint>> protonHits;
      for (int h = 0; h < evt->GetNumHits(); h++) {
         auto &hit = evt->GetHit(h);
         auto &mc = hit.GetMCSimPointArray();
         if (!mc.empty() && mc[0].trackID == 1) {
            protonHits.push_back({hit.GetPosition().Z(), hit.GetPosition()});
         }
      }
      // Sort by Z descending (vertex = high Z in digi frame)
      std::sort(protonHits.begin(), protonHits.end(),
                [](auto &a, auto &b) { return a.first > b.first; });

      if (protonHits.size() >= 3) {
         // Direction from first to 5th hit, convert to lab frame
         int nUse = std::min(5, (int)protonHits.size());
         auto p0 = protonHits[0].second;
         auto p1 = protonHits[nUse - 1].second;
         double dx = p1.X() - p0.X();
         double dy = p1.Y() - p0.Y();
         double dz = (1000.0 - p1.Z()) - (1000.0 - p0.Z()); // lab frame
         double rxy = std::sqrt(dx * dx + dy * dy);
         truthTheta = std::atan2(rxy, -dz) * 180.0 / M_PI; // lab frame angle
      }

      if (truthTheta < 0) continue;
      nTruth++;

      // --- Match Smooth3D fitted track by closest theta ---
      tSmooth->GetEntry(iEv);
      if (sArr && sArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)sArr->At(0);
         auto &ft = te->GetFittedTracks();

         double bestDth = 999;
         int bestIdx = -1;
         for (int j = 0; j < (int)ft.size(); j++) {
            auto kin = ft[j]->GetKinematics();
            if (kin.kineticEnergy <= 0) continue;
            double fitTh = kin.theta * 180.0 / M_PI;
            double dth = std::abs(fitTh - truthTheta);
            if (dth < bestDth) { bestDth = dth; bestIdx = j; }
         }

         if (bestIdx >= 0 && bestDth < 30.0) {
            auto kin = ft[bestIdx]->GetKinematics();
            double fitKE = kin.kineticEnergy;
            double fitTh = kin.theta * 180.0 / M_PI;
            double dKE_pct = (fitKE - truthKE) / truthKE * 100.0;
            double dTh = fitTh - truthTheta;

            hS_dKE->Fill(dKE_pct);
            hS_dTh->Fill(dTh);
            hS_KEvRes->Fill(truthKE, dKE_pct);
            hS_KE_vs_truth->Fill(truthKE, fitKE);
            nS++;
         }
      }

      // --- Match GNN fitted track by closest theta ---
      tGNN->GetEntry(iEv);
      if (gArr && gArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)gArr->At(0);
         auto &ft = te->GetFittedTracks();

         double bestDth = 999;
         int bestIdx = -1;
         for (int j = 0; j < (int)ft.size(); j++) {
            auto kin = ft[j]->GetKinematics();
            if (kin.kineticEnergy <= 0) continue;
            double fitTh = kin.theta * 180.0 / M_PI;
            double dth = std::abs(fitTh - truthTheta);
            if (dth < bestDth) { bestDth = dth; bestIdx = j; }
         }

         if (bestIdx >= 0 && bestDth < 30.0) {
            auto kin = ft[bestIdx]->GetKinematics();
            double fitKE = kin.kineticEnergy;
            double fitTh = kin.theta * 180.0 / M_PI;
            double dKE_pct = (fitKE - truthKE) / truthKE * 100.0;
            double dTh = fitTh - truthTheta;

            hG_dKE->Fill(dKE_pct);
            hG_dTh->Fill(dTh);
            hG_KEvRes->Fill(truthKE, dKE_pct);
            hG_KE_vs_truth->Fill(truthKE, fitKE);
            nG++;
         }
      }
   }

   // ─── Plot ───
   auto *c = new TCanvas("c", "Resolution Comparison", 1800, 1200);
   c->Divide(3, 3);

   // Row 1: KE resolution
   c->cd(1);
   hS_dKE->SetLineColor(kBlue); hS_dKE->SetLineWidth(2);
   hG_dKE->SetLineColor(kRed); hG_dKE->SetLineWidth(2);
   hS_dKE->Draw(); hG_dKE->Draw("SAME");
   auto *leg1 = new TLegend(0.15, 0.7, 0.55, 0.9);
   leg1->AddEntry(hS_dKE, Form("Smooth3D (N=%d)", nS), "l");
   leg1->AddEntry(hG_dKE, Form("GNN (N=%d)", nG), "l");
   leg1->Draw();

   // Row 1: theta resolution
   c->cd(2);
   hS_dTh->SetLineColor(kBlue); hS_dTh->SetLineWidth(2);
   hG_dTh->SetLineColor(kRed); hG_dTh->SetLineWidth(2);
   hS_dTh->Draw(); hG_dTh->Draw("SAME");
   auto *leg2 = new TLegend(0.15, 0.7, 0.55, 0.9);
   leg2->AddEntry(hS_dTh, "Smooth3D", "l");
   leg2->AddEntry(hG_dTh, "GNN", "l");
   leg2->Draw();

   // Row 1: summary
   c->cd(3);
   TPaveText *pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
   pt->SetTextAlign(12); pt->SetTextSize(0.05); pt->SetFillColor(0);
   pt->AddText(Form("#bf{Resolution Summary} (%d truth events)", nTruth));
   pt->AddText("");
   pt->AddText("#bf{Smooth3D}");
   pt->AddText(Form("  Matched: %d (%.1f%%)", nS, 100.0 * nS / nTruth));
   pt->AddText(Form("  KE: %.1f #pm %.1f %%", hS_dKE->GetMean(), hS_dKE->GetRMS()));
   pt->AddText(Form("  #theta: %.2f #pm %.2f deg", hS_dTh->GetMean(), hS_dTh->GetRMS()));
   pt->AddText("");
   pt->AddText("#bf{GNN + Smooth3D}");
   pt->AddText(Form("  Matched: %d (%.1f%%)", nG, 100.0 * nG / nTruth));
   pt->AddText(Form("  KE: %.1f #pm %.1f %%", hG_dKE->GetMean(), hG_dKE->GetRMS()));
   pt->AddText(Form("  #theta: %.2f #pm %.2f deg", hG_dTh->GetMean(), hG_dTh->GetRMS()));
   pt->Draw();

   // Row 2: KE fit vs truth
   c->cd(4);
   hS_KE_vs_truth->Draw("COLZ");
   TLine *l1 = new TLine(0, 0, 15, 15); l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->Draw();

   c->cd(5);
   hG_KE_vs_truth->Draw("COLZ");
   TLine *l2 = new TLine(0, 0, 15, 15); l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->Draw();

   // Row 2: overlay KE vs truth as profiles
   c->cd(6);
   auto *pS = hS_KE_vs_truth->ProfileX("pS"); pS->SetLineColor(kBlue); pS->SetLineWidth(2);
   auto *pG = hG_KE_vs_truth->ProfileX("pG"); pG->SetLineColor(kRed); pG->SetLineWidth(2);
   pS->SetTitle("KE profile;KE_{truth} [MeV];#LT KE_{fit} #GT [MeV]");
   pS->Draw(); pG->Draw("SAME");
   TLine *l3 = new TLine(0, 0, 15, 15); l3->SetLineColor(kBlack); l3->SetLineStyle(2); l3->Draw();
   auto *leg3 = new TLegend(0.15, 0.7, 0.5, 0.9);
   leg3->AddEntry(pS, "Smooth3D", "l");
   leg3->AddEntry(pG, "GNN", "l");
   leg3->Draw();

   // Row 3: resolution vs truth KE
   c->cd(7);
   hS_KEvRes->Draw("COLZ");

   c->cd(8);
   hG_KEvRes->Draw("COLZ");

   // Row 3: resolution profiles
   c->cd(9);
   auto *prS = hS_KEvRes->ProfileX("prS"); prS->SetLineColor(kBlue); prS->SetLineWidth(2);
   auto *prG = hG_KEvRes->ProfileX("prG"); prG->SetLineColor(kRed); prG->SetLineWidth(2);
   prS->SetTitle("KE bias vs truth;KE_{truth} [MeV];#LT #DeltaKE/KE #GT [%]");
   prS->SetMinimum(-30); prS->SetMaximum(30);
   prS->Draw(); prG->Draw("SAME");
   TLine *l4 = new TLine(0, 0, 15, 0); l4->SetLineColor(kBlack); l4->SetLineStyle(2); l4->Draw();
   auto *leg4 = new TLegend(0.55, 0.7, 0.9, 0.9);
   leg4->AddEntry(prS, "Smooth3D", "l");
   leg4->AddEntry(prG, "GNN", "l");
   leg4->Draw();

   c->SaveAs("data/gnn_smooth/resolution_comparison.png");
   std::cout << "\nSaved: data/gnn_smooth/resolution_comparison.png" << std::endl;

   fDigi->Close(); fSmooth->Close(); fGNN->Close();
}
