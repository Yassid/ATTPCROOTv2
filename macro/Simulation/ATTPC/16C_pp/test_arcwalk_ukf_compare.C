/// @file test_arcwalk_ukf_compare.C
/// @brief Compare UKF KE/theta resolution for Smooth3D vs Arc-walk.
///
/// Run after:
///   root -b -q 'test_arcwalk_ukf.C(500, false)'
///   root -b -q 'test_arcwalk_ukf.C(500, true)'
///
/// Then: root -b -q test_arcwalk_ukf_compare.C

void test_arcwalk_ukf_compare()
{
   TString inOutDir = "./data/";

   auto fMC = TFile::Open(inOutDir + "attpcsim.root");
   auto fSmooth = TFile::Open(inOutDir + "output_smooth3d_ukf.root");
   auto fArcWalk = TFile::Open(inOutDir + "output_arcwalk_ukf.root");

   if (!fMC || !fSmooth || !fArcWalk) {
      std::cerr << "Run test_arcwalk_ukf.C with both modes first." << std::endl;
      return;
   }

   auto tMC = (TTree *)fMC->Get("cbmsim");
   auto tSmooth = (TTree *)fSmooth->Get("cbmsim");
   auto tArcWalk = (TTree *)fArcWalk->Get("cbmsim");

   TClonesArray *mcTracks = nullptr;
   tMC->SetBranchAddress("MCTrack", &mcTracks);

   TClonesArray *trkSmooth = nullptr, *trkArcWalk = nullptr;
   tSmooth->SetBranchAddress("AtTrackingEvent", &trkSmooth);
   tArcWalk->SetBranchAddress("AtTrackingEvent", &trkArcWalk);

   int nEntries = std::min({(int)tMC->GetEntries(), (int)tSmooth->GetEntries(), (int)tArcWalk->GetEntries()});

   // KE error histograms
   TH1F *hKE_s = new TH1F("hKE_s", "KE error (Smooth3D);(KE_{reco} - KE_{true}) / KE_{true} [%];Events", 60, -30, 30);
   TH1F *hKE_a = new TH1F("hKE_a", "KE error (Arc-walk);(KE_{reco} - KE_{true}) / KE_{true} [%];Events", 60, -30, 30);
   TH1F *hTh_s =
      new TH1F("hTh_s", "#theta error (Smooth3D);#theta_{reco} - #theta_{true} [deg];Events", 60, -10, 10);
   TH1F *hTh_a =
      new TH1F("hTh_a", "#theta error (Arc-walk);#theta_{reco} - #theta_{true} [deg];Events", 60, -10, 10);
   TH1F *hNCl_s = new TH1F("hNCl_s", "Clusters (Smooth3D);N clusters;Events", 60, 0, 60);
   TH1F *hNCl_a = new TH1F("hNCl_a", "Clusters (Arc-walk);N clusters;Events", 60, 0, 60);

   int nFit_s = 0, nFit_a = 0;

   for (int i = 0; i < nEntries; i++) {
      tMC->GetEntry(i);
      tSmooth->GetEntry(i);
      tArcWalk->GetEntry(i);

      // Find proton MC truth
      double keTrue = -1, thetaTrue = -1;
      for (int j = 0; j < mcTracks->GetEntries(); j++) {
         auto *t = (AtMCTrack *)mcTracks->At(j);
         if (t->GetPdgCode() == 2212) {
            double px = t->GetPx() * 1e3, py = t->GetPy() * 1e3, pz = t->GetPz() * 1e3;
            double p = std::sqrt(px * px + py * py + pz * pz);
            keTrue = std::sqrt(p * p + 938.272 * 938.272) - 938.272;
            thetaTrue = std::acos(pz / p) * 180.0 / M_PI;
            break;
         }
      }
      if (keTrue < 0)
         continue;

      // Extract fitted KE from Smooth3D
      if (trkSmooth && trkSmooth->GetEntries() > 0) {
         auto *evt = (AtTrackingEvent *)trkSmooth->At(0);
         auto &fitted = evt->GetFittedTracks();
         if (!fitted.empty()) {
            try {
               auto &kin = fitted[0]->GetKinematics();
               if (kin.kineticEnergy > 0) {
                  nFit_s++;
                  hKE_s->Fill((kin.kineticEnergy - keTrue) / keTrue * 100);
                  hTh_s->Fill(kin.theta * 180.0 / M_PI - thetaTrue);
               }
            } catch (...) {
            }
         }
      }

      // Extract fitted KE from Arc-walk
      if (trkArcWalk && trkArcWalk->GetEntries() > 0) {
         auto *evt = (AtTrackingEvent *)trkArcWalk->At(0);
         auto &fitted = evt->GetFittedTracks();
         if (!fitted.empty()) {
            try {
               auto &kin = fitted[0]->GetKinematics();
               if (kin.kineticEnergy > 0) {
                  nFit_a++;
                  hKE_a->Fill((kin.kineticEnergy - keTrue) / keTrue * 100);
                  hTh_a->Fill(kin.theta * 180.0 / M_PI - thetaTrue);
               }
            } catch (...) {
            }
         }
      }
   }

   std::cout << "\n========================================" << std::endl;
   std::cout << "   Smooth3D vs Arc-Walk — UKF Resolution" << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "Events:              " << nEntries << std::endl;
   printf("                     %-12s%-12s\n", "Smooth3D", "Arc-Walk");
   printf("Fitted tracks:       %-12d%-12d\n", nFit_s, nFit_a);
   printf("KE bias [%%]:         %-12.2f%-12.2f\n", hKE_s->GetMean(), hKE_a->GetMean());
   printf("KE RMS [%%]:          %-12.2f%-12.2f\n", hKE_s->GetRMS(), hKE_a->GetRMS());
   printf("Theta bias [deg]:    %-12.3f%-12.3f\n", hTh_s->GetMean(), hTh_a->GetMean());
   printf("Theta RMS [deg]:     %-12.3f%-12.3f\n", hTh_s->GetRMS(), hTh_a->GetRMS());
   std::cout << "========================================" << std::endl;

   // Plot
   TCanvas *c = new TCanvas("cUKFComp", "UKF Resolution Comparison", 1200, 800);
   c->Divide(2, 2);

   c->cd(1);
   hKE_s->SetLineColor(kBlue);
   hKE_s->SetLineWidth(2);
   hKE_a->SetLineColor(kRed);
   hKE_a->SetLineWidth(2);
   hKE_s->Draw("hist");
   hKE_a->Draw("hist same");
   auto l1 = new TLegend(0.15, 0.7, 0.48, 0.88);
   l1->AddEntry(hKE_s, Form("Smooth3D (RMS %.1f%%)", hKE_s->GetRMS()), "l");
   l1->AddEntry(hKE_a, Form("Arc-walk (RMS %.1f%%)", hKE_a->GetRMS()), "l");
   l1->Draw();

   c->cd(2);
   hTh_s->SetLineColor(kBlue);
   hTh_s->SetLineWidth(2);
   hTh_a->SetLineColor(kRed);
   hTh_a->SetLineWidth(2);
   hTh_s->Draw("hist");
   hTh_a->Draw("hist same");
   auto l2 = new TLegend(0.15, 0.7, 0.48, 0.88);
   l2->AddEntry(hTh_s, Form("Smooth3D (RMS %.2f#circ)", hTh_s->GetRMS()), "l");
   l2->AddEntry(hTh_a, Form("Arc-walk (RMS %.2f#circ)", hTh_a->GetRMS()), "l");
   l2->Draw();

   // KE error vs true KE scatter
   c->cd(3);
   // Placeholder: KE distributions
   TH1F *hKETrue = new TH1F("hKETrue", "KE_{true};KE [MeV];Events", 50, 0, 50);
   for (int i = 0; i < nEntries; i++) {
      tMC->GetEntry(i);
      for (int j = 0; j < mcTracks->GetEntries(); j++) {
         auto *t = (AtMCTrack *)mcTracks->At(j);
         if (t->GetPdgCode() == 2212) {
            double px = t->GetPx() * 1e3, py = t->GetPy() * 1e3, pz = t->GetPz() * 1e3;
            double p = std::sqrt(px * px + py * py + pz * pz);
            hKETrue->Fill(std::sqrt(p * p + 938.272 * 938.272) - 938.272);
            break;
         }
      }
   }
   hKETrue->Draw("hist");

   c->cd(4);
   // Efficiency text
   TPaveText *pt = new TPaveText(0.1, 0.3, 0.9, 0.7, "NDC");
   pt->AddText(Form("Smooth3D: %d/%d fitted (%.1f%%)", nFit_s, nEntries / 2, 100.0 * nFit_s / (nEntries / 2)));
   pt->AddText(Form("Arc-walk: %d/%d fitted (%.1f%%)", nFit_a, nEntries / 2, 100.0 * nFit_a / (nEntries / 2)));
   pt->AddText(Form("KE RMS: %.1f%% vs %.1f%%", hKE_s->GetRMS(), hKE_a->GetRMS()));
   pt->AddText(Form("#theta RMS: %.2f#circ vs %.2f#circ", hTh_s->GetRMS(), hTh_a->GetRMS()));
   pt->Draw();

   c->SaveAs("data/arcwalk_ukf_comparison.png");
   std::cout << "Saved: data/arcwalk_ukf_comparison.png" << std::endl;

   fMC->Close();
   fSmooth->Close();
   fArcWalk->Close();
}
