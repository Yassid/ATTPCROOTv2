/// @file show_kinematics.C
/// @brief Display kinematics results interactively.
/// Run: root -l show_kinematics.C

void show_kinematics()
{
   TFile *fm = TFile::Open("data/attpcsim.root");
   TFile *ff = TFile::Open("data/output_ukf_only.root");
   if (!fm || !ff) {
      std::cout << "Missing files." << std::endl;
      return;
   }

   TTree *tm = (TTree *)fm->Get("cbmsim");
   TTree *tf = (TTree *)ff->Get("cbmsim");
   TClonesArray *ma = nullptr, *fa = nullptr;
   tm->SetBranchAddress("MCTrack", &ma);
   tf->SetBranchAddress("AtTrackingEvent", &fa);

   int nEvents = std::min(tm->GetEntries(), tf->GetEntries());

   std::vector<double> truthTheta, truthKE;
   std::vector<double> recoTheta, recoKE;
   std::vector<double> keErr, thErr;
   std::vector<double> vtxX, vtxY;       // reconstructed vertex XY
   std::vector<int> recoEventId;         // event number for each reco point

   for (int i = 0; i < nEvents; i++) {
      tm->GetEntry(i);
      tf->GetEntry(i);

      double pT = -1, thT = -1, keT = -1;
      for (int j = 0; j < ma->GetEntries(); j++) {
         AtMCTrack *m = (AtMCTrack *)ma->At(j);
         if (m->GetPdgCode() == 2212) {
            ROOT::Math::XYZVector mom(m->GetPx(), m->GetPy(), m->GetPz());
            pT = mom.R() * 1e3;
            thT = mom.Theta() * 180.0 / M_PI;
            keT = std::sqrt(pT * pT + 938.272 * 938.272) - 938.272;
            break;
         }
      }
      if (pT < 0) continue;

      truthTheta.push_back(thT);
      truthKE.push_back(keT);

      AtTrackingEvent *te = (AtTrackingEvent *)fa->At(0);
      auto &ft = te->GetFittedTracks();
      if (ft.empty()) continue;

      double keR = ft[0]->GetKinematics().kineticEnergy;
      double thR = ft[0]->GetKinematics().theta * 180.0 / M_PI;
      if (keR < 0.05 || keR > 50 || thR < 10 || thR > 170) continue;

      recoTheta.push_back(thR);
      recoKE.push_back(keR);
      keErr.push_back((keR - keT) / keT * 100);
      thErr.push_back(thR - thT);

      auto vtx = ft[0]->GetVertex();
      vtxX.push_back(vtx.X());
      vtxY.push_back(vtx.Y());
      recoEventId.push_back(i);
   }

   // Canvas 1: Overlay
   TCanvas *c1 = new TCanvas("c1", "Kinematic Curve", 800, 600);
   TGraph *gT = new TGraph(truthTheta.size(), truthTheta.data(), truthKE.data());
   gT->SetTitle("16C(p,p): Truth vs Reconstructed;#theta_{lab} [deg];KE [MeV]");
   gT->SetMarkerStyle(20);
   gT->SetMarkerSize(0.4);
   gT->SetMarkerColor(kBlack);
   gT->Draw("AP");

   TGraph *gR = new TGraph(recoTheta.size(), recoTheta.data(), recoKE.data());
   gR->SetMarkerStyle(20);
   gR->SetMarkerSize(0.4);
   gR->SetMarkerColor(kRed);
   gR->Draw("P SAME");

   TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg->AddEntry(gT, Form("MC Truth (%d)", (int)truthTheta.size()), "p");
   leg->AddEntry(gR, Form("UKF Reco (%d)", (int)recoTheta.size()), "p");
   leg->Draw();

   // Canvas 2: Errors
   TCanvas *c2 = new TCanvas("c2", "Reconstruction Errors", 1000, 500);
   c2->Divide(2, 1);

   c2->cd(1);
   TH1F *hKE = new TH1F("hKE", "KE error;(KE_{reco} - KE_{true}) / KE_{true} [%];Events", 80, -50, 50);
   for (double e : keErr) hKE->Fill(e);
   hKE->Fit("gaus", "Q");
   hKE->Draw("hist");
   if (hKE->GetFunction("gaus")) hKE->GetFunction("gaus")->Draw("same");

   c2->cd(2);
   TH1F *hTh = new TH1F("hTh", "#theta error;#theta_{reco} - #theta_{true} [deg];Events", 80, -30, 30);
   for (double e : thErr) hTh->Fill(e);
   hTh->Fit("gaus", "Q");
   hTh->Draw("hist");
   if (hTh->GetFunction("gaus")) hTh->GetFunction("gaus")->Draw("same");

   // Canvas 3: Low energy zoom with event labels
   TCanvas *c3 = new TCanvas("c3", "Low Energy Detail", 800, 600);
   // Plot only events with KE_true < 5 MeV
   std::vector<double> lowThT, lowKET, lowThR, lowKER;
   std::vector<int> lowEvId;
   for (size_t k = 0; k < recoTheta.size(); k++) {
      // Find the matching truth for this reco event
      int evId = recoEventId[k];
      // Find truth index
      double thT = -1, keT = -1;
      for (size_t m = 0; m < truthTheta.size(); m++) {
         // truthTheta has all proton events, recoEventId maps to absolute event number
         // We need to find the truth that corresponds to this event
         // Since truth is filled for all proton events in order, and reco is a subset,
         // we need the actual truth values at this event
         break;
      }
   }
   // Simpler: just re-read for low energy events
   {
      tm->SetBranchAddress("MCTrack", &ma);
      tf->SetBranchAddress("AtTrackingEvent", &fa);
      for (size_t k = 0; k < recoEventId.size(); k++) {
         int evId = recoEventId[k];
         tm->GetEntry(evId);
         tf->GetEntry(evId);
         double pT = -1, thT = -1, keT = -1;
         for (int j = 0; j < ma->GetEntries(); j++) {
            AtMCTrack *m = (AtMCTrack *)ma->At(j);
            if (m->GetPdgCode() == 2212) {
               ROOT::Math::XYZVector mom(m->GetPx(), m->GetPy(), m->GetPz());
               pT = mom.R() * 1e3;
               thT = mom.Theta() * 180.0 / M_PI;
               keT = std::sqrt(pT * pT + 938.272 * 938.272) - 938.272;
               break;
            }
         }
         if (keT < 5.0 && keT > 0) {
            lowThT.push_back(thT);
            lowKET.push_back(keT);
            lowThR.push_back(recoTheta[k]);
            lowKER.push_back(recoKE[k]);
            lowEvId.push_back(evId);
         }
      }
   }
   TGraph *gLowT = new TGraph(lowThT.size(), lowThT.data(), lowKET.data());
   gLowT->SetTitle("Low Energy (<5 MeV): Truth (black) vs Reco (red);#theta_{lab} [deg];KE [MeV]");
   gLowT->SetMarkerStyle(20);
   gLowT->SetMarkerSize(0.6);
   gLowT->SetMarkerColor(kBlack);
   gLowT->Draw("AP");
   TGraph *gLowR = new TGraph(lowThR.size(), lowThR.data(), lowKER.data());
   gLowR->SetMarkerStyle(20);
   gLowR->SetMarkerSize(0.6);
   gLowR->SetMarkerColor(kRed);
   gLowR->Draw("P SAME");
   // Draw lines connecting truth to reco for deviating events
   for (size_t k = 0; k < lowEvId.size(); k++) {
      double err = (lowKER[k] - lowKET[k]) / lowKET[k] * 100;
      if (std::abs(err) > 15) {
         TLine *line = new TLine(lowThT[k], lowKET[k], lowThR[k], lowKER[k]);
         line->SetLineColor(kMagenta);
         line->Draw();
         TLatex *label = new TLatex(lowThR[k] + 0.5, lowKER[k], Form("%d", lowEvId[k]));
         label->SetTextSize(0.025);
         label->SetTextColor(kMagenta);
         label->Draw();
      }
   }

   // Print deviating low-energy events to terminal
   std::cout << "\n=== Deviating low-energy events (KE<5 MeV, |err|>15%) ===" << std::endl;
   std::cout << "Event  KE_true  KE_reco  err%     theta_true  theta_reco  vtxR" << std::endl;
   for (size_t k = 0; k < lowEvId.size(); k++) {
      double err = (lowKER[k] - lowKET[k]) / lowKET[k] * 100;
      if (std::abs(err) > 15) {
         double vR = std::sqrt(vtxX[k] * vtxX[k] + vtxY[k] * vtxY[k]);
         std::cout << lowEvId[k] << "    " << lowKET[k] << "    " << lowKER[k] << "    " << err << "    "
                   << lowThT[k] << "    " << lowThR[k] << "    " << vR << std::endl;
      }
   }

   // Canvas 4: Vertex XY
   TCanvas *c4 = new TCanvas("c4", "Vertex XY", 600, 600);
   TGraph *gVtx = new TGraph(vtxX.size(), vtxX.data(), vtxY.data());
   gVtx->SetTitle("Reconstructed Vertex;X [mm];Y [mm]");
   gVtx->SetMarkerStyle(20);
   gVtx->SetMarkerSize(0.4);
   gVtx->SetMarkerColor(kBlue);
   gVtx->Draw("AP");
   // Draw beam axis marker
   TMarker *beamMark = new TMarker(0, 0, 29);
   beamMark->SetMarkerSize(2);
   beamMark->SetMarkerColor(kRed);
   beamMark->Draw();

   std::cout << "Proton events: " << truthTheta.size() << ", Reconstructed: " << recoTheta.size() << std::endl;
}
