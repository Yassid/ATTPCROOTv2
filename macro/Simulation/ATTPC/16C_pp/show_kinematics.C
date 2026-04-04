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
   std::vector<double> vtxX, vtxY; // reconstructed vertex XY

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

   // Canvas 3: Vertex XY
   TCanvas *c3 = new TCanvas("c3", "Vertex Position", 600, 600);
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
