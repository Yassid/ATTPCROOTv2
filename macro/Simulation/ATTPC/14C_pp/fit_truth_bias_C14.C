/// @file fit_truth_bias_C14.C
/// @brief Fitted-vs-true KE and theta for the recoil proton, UKF and GENFIT, from the sim.
///
/// The data shows the reconstructed proton KE sitting +8 to +14 % above elastic two-body
/// kinematics for KE 14-25 MeV, in BOTH fitters, which is what pushes the E_x locus to -1.4 MeV
/// at backward theta_cm and destroys the elastic yield there. Two families of cause:
///   * a calibration error in the DATA (B field, drift velocity, vertex, energy-loss model), in
///     which case the simulation -- digitised and reconstructed with consistent parameters --
///     will show no such bias;
///   * an algorithmic bias in the fitters, which the simulation WILL reproduce.
/// This macro measures it in the sim with MC truth, so the two can be told apart.
///
/// Truth matching is the same as acceptance_C14.C: the primary proton from MCTrack, paired with
/// the fitted track in the same entry that is closest in theta, requiring |dtheta| < 10 deg so a
/// mismatched track cannot masquerade as a bias.
///
///   root -b -q 'fit_truth_bias_C14.C("/mnt/f/a1954_C14_acc/","ukf",kFALSE)'
///   root -b -q 'fit_truth_bias_C14.C("/mnt/f/a1954_C14_acc_gf/","genfit",kTRUE)'

void fit_truth_bias_C14(TString accDir = "/mnt/f/a1954_C14_acc/", TString fitter = "ukf", Bool_t useXtr = kFALSE,
                        TString tag = "gs", Int_t s0 = 1001, Int_t s1 = 1005, Double_t chi2Cut = 1e9)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, m_p = 1.007825 * u;

   // KE bias in bins of TRUE KE, and theta bias in bins of true theta
   const int NK = 16;
   const double kLo = 0, kHi = 40;
   auto *pKE = new TProfile("pKE", "", NK, kLo, kHi, "s");
   auto *pRel = new TProfile("pRel", "", NK, kLo, kHi, "s");
   auto *pTh = new TProfile("pTh", "", 18, 0, 90, "s");
   auto *hKE2 = new TH2D("hKE2", "", 80, kLo, kHi, 80, -6, 6);
   auto *hTh2 = new TH2D("hTh2", "", 90, 0, 90, 80, -10, 10);
   long nPair = 0;

   for (int s = s0; s <= s1; ++s) {
      TString sf = TString::Format("%s%s_s%d_sim.root", accDir.Data(), tag.Data(), s);
      TString ff = TString::Format("%s%s_s%d_%s.root", accDir.Data(), tag.Data(), s, fitter.Data());
      TFile *fs = TFile::Open(sf), *ff2 = TFile::Open(ff);
      if (!fs || fs->IsZombie() || !ff2 || ff2->IsZombie()) {
         printf("\033[1;31mmissing seed %d (%s / %s)\033[0m\n", s, sf.Data(), ff.Data());
         continue;
      }
      TTree *ts = (TTree *)fs->Get("cbmsim"), *tf = (TTree *)ff2->Get("cbmsim");
      if (!ts || !tf || ts->GetEntries() != tf->GetEntries()) {
         printf("\033[1;31mseed %d: entry mismatch\033[0m\n", s);
         continue;
      }
      TClonesArray *mc = nullptr, *te = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);
      Long64_t N = ts->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *tr = (AtMCTrack *)mc->At(k);
            if (!tr || tr->GetPdgCode() != 2212 || tr->GetMotherId() != -1)
               continue;
            double px = tr->GetPx() * 1000, py = tr->GetPy() * 1000, pz = tr->GetPz() * 1000;
            double p = std::sqrt(px * px + py * py + pz * pz);
            if (p <= 0)
               continue;
            keT = std::sqrt(p * p + m_p * m_p) - m_p;
            thT = std::acos(pz / p) * TMath::RadToDeg();
            break;
         }
         if (keT <= 0)
            continue;
         tf->GetEntry(i);
         if (!te || te->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         // closest-in-theta fitted track, then the 10 deg truth-match gate
         double bestD = 1e9, bKE = -1, bTh = -1;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            const auto &md = ft->GetTrackMetadata();
            double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
            if (!(ndf > 0 && chi2 / ndf < chi2Cut))
               continue;
            const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
            double ke = kin.kineticEnergy, th = kin.theta * TMath::RadToDeg();
            if (!(ke > 0 && ke < 1000))
               continue;
            double d = std::fabs(th - thT);
            if (d < bestD) {
               bestD = d;
               bKE = ke;
               bTh = th;
            }
         }
         if (bKE < 0 || bestD > 10.0)
            continue;
         ++nPair;
         pKE->Fill(keT, bKE - keT);
         pRel->Fill(keT, 100.0 * (bKE - keT) / keT);
         pTh->Fill(thT, bTh - thT);
         hKE2->Fill(keT, bKE - keT);
         hTh2->Fill(thT, bTh - thT);
      }
      fs->Close();
      ff2->Close();
   }

   printf("\n===== %s on the simulation: %ld truth-matched protons =====\n", fitter.Data(), nPair);
   printf("  KE_true [MeV] |   <KE_fit-KE_true>   rms   |  rel bias [%%]   N\n");
   for (int b = 1; b <= NK; ++b) {
      if (pKE->GetBinEntries(b) < 20)
         continue;
      printf("  %5.1f-%5.1f   |  %+10.3f  %8.3f  |  %+8.2f     %6.0f\n", pKE->GetBinLowEdge(b),
             pKE->GetBinLowEdge(b) + pKE->GetBinWidth(b), pKE->GetBinContent(b), pKE->GetBinError(b),
             pRel->GetBinContent(b), pKE->GetBinEntries(b));
   }
   printf("\n  theta_true [deg] | <theta_fit-theta_true>   rms       N\n");
   for (int b = 1; b <= pTh->GetNbinsX(); ++b) {
      if (pTh->GetBinEntries(b) < 20)
         continue;
      printf("  %5.0f-%5.0f      |  %+12.3f  %10.3f  %7.0f\n", pTh->GetBinLowEdge(b),
             pTh->GetBinLowEdge(b) + pTh->GetBinWidth(b), pTh->GetBinContent(b), pTh->GetBinError(b),
             pTh->GetBinEntries(b));
   }

   TCanvas *c1 = new TCanvas("c1", "truth bias", 1500, 950);
   c1->Divide(2, 2);
   c1->cd(1);
   hKE2->SetTitle(TString::Format("%s: KE_{fit} - KE_{true};KE_{true} [MeV];#DeltaKE [MeV]", fitter.Data()));
   hKE2->Draw("colz");
   pKE->SetLineColor(kRed + 1);
   pKE->SetLineWidth(3);
   pKE->SetMarkerStyle(20);
   pKE->SetMarkerColor(kRed + 1);
   pKE->Draw("same");
   c1->cd(2);
   pRel->SetTitle(TString::Format("%s: relative KE bias;KE_{true} [MeV];bias [%%]", fitter.Data()));
   pRel->SetLineColor(kAzure + 2);
   pRel->SetLineWidth(3);
   pRel->SetMarkerStyle(20);
   pRel->GetYaxis()->SetRangeUser(-20, 20);
   pRel->Draw();
   auto *z = new TLine(kLo, 0, kHi, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   c1->cd(3);
   hTh2->SetTitle(TString::Format("%s: #theta_{fit} - #theta_{true};#theta_{true} [deg];#Delta#theta [deg]",
                                  fitter.Data()));
   hTh2->Draw("colz");
   pTh->SetLineColor(kRed + 1);
   pTh->SetLineWidth(3);
   pTh->SetMarkerStyle(20);
   pTh->SetMarkerColor(kRed + 1);
   pTh->Draw("same");
   c1->cd(4);
   pTh->SetTitle(TString::Format("%s: #theta bias;#theta_{true} [deg];#Delta#theta [deg]", fitter.Data()));
   pTh->GetYaxis()->SetRangeUser(-4, 4);
   pTh->Draw();
   auto *z2 = new TLine(0, 0, 90, 0);
   z2->SetLineStyle(2);
   z2->SetLineColor(kGray + 2);
   z2->Draw();

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString png = here + "/diagnostics/fit_truth_bias_" + fitter + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/diagnostics/fit_truth_bias_" + fitter + ".root", "RECREATE");
   pKE->Write("dKE");
   pRel->Write("relKE");
   pTh->Write("dTheta");
   hKE2->Write("dKE2D");
   hTh2->Write("dTheta2D");
   fo.Close();
   printf("\nwrote %s\n\n", png.Data());
}
