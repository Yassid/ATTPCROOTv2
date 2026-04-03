/// @file analyze_kinematics.C
/// @brief Analyze angle-energy kinematics: MC truth vs UKF reconstructed.
///
/// Reads MC truth (attpcsim.root), digitized (output_digi.root), and
/// fitted (output_ukf_only.root) data. Plots the kinematic curve
/// (theta_lab vs KE) for both truth and reconstructed protons.
///
/// Run: root -b -q analyze_kinematics.C
///      root -l analyze_kinematics.C   // interactive

void analyze_kinematics()
{
   FairLogger::GetLogger()->SetLogScreenLevel("FATAL");

   TFile *fMC = TFile::Open("data/attpcsim.root");
   TFile *fDigi = TFile::Open("data/output_digi.root");
   TFile *fFit = TFile::Open("data/output_ukf_only.root");

   if (!fMC || !fDigi || !fFit) {
      std::cerr << "Missing files. Run sim, digi, and ukf_only first." << std::endl;
      return;
   }

   TTree *tMC = (TTree *)fMC->Get("cbmsim");
   TTree *tDigi = (TTree *)fDigi->Get("cbmsim");
   TTree *tFit = (TTree *)fFit->Get("cbmsim");

   TClonesArray *mcTracks = nullptr;
   TClonesArray *patEvtArr = nullptr;
   TClonesArray *fitEvtArr = nullptr;

   tMC->SetBranchAddress("MCTrack", &mcTracks);
   tDigi->SetBranchAddress("AtPatternEvent", &patEvtArr);
   tFit->SetBranchAddress("AtTrackingEvent", &fitEvtArr);

   int nEvents = std::min({(int)tMC->GetEntries(), (int)tDigi->GetEntries(), (int)tFit->GetEntries()});

   // Histograms
   TH2F *hTruth = new TH2F("hTruth", "MC Truth;#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 100, 0, 50);
   TH2F *hReco = new TH2F("hReco", "Reconstructed (UKF);#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 100, 0, 50);
   TH2F *hBrho = new TH2F("hBrho", "Brho seed;#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 100, 0, 50);

   // Graphs for scatter plot
   std::vector<double> truthTheta, truthKE;
   std::vector<double> recoTheta, recoKE;
   std::vector<double> brhoTheta, brhoKE;

   // 1D distributions
   TH1F *hKEerr = new TH1F("hKEerr", "KE error;(KE_{reco} - KE_{true}) / KE_{true} [%];Events", 100, -50, 50);
   TH1F *hThetaErr = new TH1F("hThetaErr", "#theta error;#theta_{reco} - #theta_{true} [deg];Events", 100, -30, 30);

   int nFitted = 0, nGood = 0;

   for (int i = 0; i < nEvents; i++) {
      tMC->GetEntry(i);
      tDigi->GetEntry(i);
      tFit->GetEntry(i);

      // MC truth proton
      double pTrue = -1, thetaTrue = -1, keTrue = -1;
      for (int j = 0; j < mcTracks->GetEntries(); j++) {
         AtMCTrack *mc = (AtMCTrack *)mcTracks->At(j);
         if (mc->GetPdgCode() == 2212) {
            ROOT::Math::XYZVector mom(mc->GetPx(), mc->GetPy(), mc->GetPz());
            pTrue = mom.R() * 1e3; // GeV → MeV
            thetaTrue = mom.Theta() * 180.0 / M_PI;
            keTrue = std::sqrt(pTrue * pTrue + 938.272 * 938.272) - 938.272;
            break;
         }
      }
      if (pTrue < 0)
         continue;

      hTruth->Fill(thetaTrue, keTrue);
      truthTheta.push_back(thetaTrue);
      truthKE.push_back(keTrue);

      // Brho seed from pattern recognition
      if (patEvtArr && patEvtArr->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)patEvtArr->At(0);
         if (!pe->GetTrackCand().empty()) {
            auto &tr = pe->GetTrackCand()[0];
            double R = tr.GetGeoRadius() / 1000.0;
            double th = tr.GetGeoTheta();
            double sinTh = std::sin(th);
            if (std::abs(sinTh) > 1e-9) {
               double pBrho = 0.3 * 2.85 * R / sinTh * 1000;
               double keBrho = std::sqrt(pBrho * pBrho + 938.272 * 938.272) - 938.272;
               // GeoTheta is in digi frame → lab: 180 - theta_digi
               double thetaBrho = 180.0 - th * 180.0 / M_PI;
               hBrho->Fill(thetaBrho, keBrho);
               brhoTheta.push_back(thetaBrho);
               brhoKE.push_back(keBrho);
            }
         }
      }

      // UKF reconstructed
      if (fitEvtArr && fitEvtArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)fitEvtArr->At(0);
         auto &ft = te->GetFittedTracks();
         if (!ft.empty()) {
            nFitted++;
            auto kin = ft[0]->GetKinematics();
            double keReco = kin.kineticEnergy;
            double thetaReco = kin.theta * 180.0 / M_PI;

            // Check if reasonable
            if (keReco > 0.05 && keReco < 50.0 && thetaReco > 10 && thetaReco < 170) {
               nGood++;
               hReco->Fill(thetaReco, keReco);
               recoTheta.push_back(thetaReco);
               recoKE.push_back(keReco);

               hKEerr->Fill((keReco - keTrue) / keTrue * 100);
               hThetaErr->Fill(thetaReco - thetaTrue);
            }
         }
      }
   }

   std::cout << "Events: " << nEvents << ", Fitted: " << nFitted << ", Good: " << nGood << std::endl;

   // === Plots ===
   TCanvas *c1 = new TCanvas("c1", "Kinematic Curves", 1400, 600);
   c1->Divide(3, 1);

   // MC Truth scatter
   c1->cd(1);
   TGraph *gTruth = new TGraph(truthTheta.size(), truthTheta.data(), truthKE.data());
   gTruth->SetTitle("MC Truth;#theta_{lab} [deg];KE [MeV]");
   gTruth->SetMarkerStyle(20);
   gTruth->SetMarkerSize(0.3);
   gTruth->SetMarkerColor(kBlack);
   gTruth->Draw("AP");

   // Brho scatter
   c1->cd(2);
   TGraph *gBrho = new TGraph(brhoTheta.size(), brhoTheta.data(), brhoKE.data());
   gBrho->SetTitle("Brho Seed;#theta_{lab} [deg];KE [MeV]");
   gBrho->SetMarkerStyle(20);
   gBrho->SetMarkerSize(0.3);
   gBrho->SetMarkerColor(kBlue);
   gBrho->Draw("AP");
   // Overlay truth as gray
   TGraph *gTruth2 = new TGraph(truthTheta.size(), truthTheta.data(), truthKE.data());
   gTruth2->SetMarkerStyle(20);
   gTruth2->SetMarkerSize(0.2);
   gTruth2->SetMarkerColor(kGray);
   gTruth2->Draw("P SAME");

   // Reco scatter
   c1->cd(3);
   TGraph *gReco = new TGraph(recoTheta.size(), recoTheta.data(), recoKE.data());
   gReco->SetTitle("Reconstructed (UKF);#theta_{lab} [deg];KE [MeV]");
   gReco->SetMarkerStyle(20);
   gReco->SetMarkerSize(0.3);
   gReco->SetMarkerColor(kRed);
   gReco->Draw("AP");
   // Overlay truth as gray
   TGraph *gTruth3 = new TGraph(truthTheta.size(), truthTheta.data(), truthKE.data());
   gTruth3->SetMarkerStyle(20);
   gTruth3->SetMarkerSize(0.2);
   gTruth3->SetMarkerColor(kGray);
   gTruth3->Draw("P SAME");

   c1->SaveAs("data/kinematics_curves.png");
   std::cout << "Saved: data/kinematics_curves.png" << std::endl;

   // Error distributions
   TCanvas *c2 = new TCanvas("c2", "Reconstruction Errors", 1000, 500);
   c2->Divide(2, 1);
   c2->cd(1);
   hKEerr->Fit("gaus", "Q");
   hKEerr->Draw("hist");
   if (hKEerr->GetFunction("gaus"))
      hKEerr->GetFunction("gaus")->Draw("same");
   c2->cd(2);
   hThetaErr->Fit("gaus", "Q");
   hThetaErr->Draw("hist");
   if (hThetaErr->GetFunction("gaus"))
      hThetaErr->GetFunction("gaus")->Draw("same");
   c2->SaveAs("data/kinematics_errors.png");
   std::cout << "Saved: data/kinematics_errors.png" << std::endl;

   // 2D overlay
   TCanvas *c3 = new TCanvas("c3", "Truth vs Reco Overlay", 700, 600);
   gTruth->SetTitle("Kinematic Curve: Truth (black) vs Reco (red);#theta_{lab} [deg];KE [MeV]");
   gTruth->Draw("AP");
   gReco->Draw("P SAME");
   TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg->AddEntry(gTruth, Form("MC Truth (%d)", (int)truthTheta.size()), "p");
   leg->AddEntry(gReco, Form("UKF Reco (%d)", (int)recoTheta.size()), "p");
   leg->Draw();
   c3->SaveAs("data/kinematics_overlay.png");
   std::cout << "Saved: data/kinematics_overlay.png" << std::endl;
}
