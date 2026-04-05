void chi2_cuts_kinematics()
{
   TFile *fm = TFile::Open("data/attpcsim.root");
   TFile *fd = TFile::Open("data/output_digi.root");
   TFile *ff = TFile::Open("data/output_ukf_only.root");
   TTree *tm = (TTree *)fm->Get("cbmsim");
   TTree *td = (TTree *)fd->Get("cbmsim");
   TTree *tf = (TTree *)ff->Get("cbmsim");
   TClonesArray *ma = nullptr, *pa = nullptr, *fa = nullptr;
   tm->SetBranchAddress("MCTrack", &ma);
   td->SetBranchAddress("AtPatternEvent", &pa);
   tf->SetBranchAddress("AtTrackingEvent", &fa);

   int nEv = std::min({(int)tm->GetEntries(), (int)td->GetEntries(), (int)tf->GetEntries()});

   struct EvData { double chi2, thTrue, keTrue, thReco, keReco, keErr; };
   std::vector<EvData> data;

   for (int i = 0; i < nEv; i++) {
      tm->GetEntry(i);
      td->GetEntry(i);
      tf->GetEntry(i);

      double pT = -1, thT = -1, keT = -1;
      for (int j = 0; j < ma->GetEntries(); j++) {
         auto *mc = (AtMCTrack *)ma->At(j);
         if (mc->GetPdgCode() == 2212) {
            ROOT::Math::XYZVector mom(mc->GetPx(), mc->GetPy(), mc->GetPz());
            pT = mom.R() * 1e3;
            thT = mom.Theta() * 180.0 / M_PI;
            keT = std::sqrt(pT * pT + 938.272 * 938.272) - 938.272;
            break;
         }
      }
      if (pT < 0) continue;

      if (!fa || fa->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)fa->At(0);
      auto &ft = te->GetFittedTracks();
      if (ft.empty()) continue;

      auto kin = ft[0]->GetKinematics();
      if (kin.kineticEnergy < 0.05 || kin.kineticEnergy > 50) continue;
      double thR = kin.theta * 180.0 / M_PI;
      if (thR < 10 || thR > 170) continue;

      auto [pval, chi2, bchi2, ndf, bndf, conv] = ft[0]->GetStats();
      if (ndf <= 0 || chi2 <= 0) continue;

      data.push_back({chi2, thT, keT, thR, kin.kineticEnergy,
                       (kin.kineticEnergy - keT) / keT * 100});
   }

   double cuts[] = {0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
   int nCuts = 6;
   int colors[] = {kBlue, kGreen+2, kOrange+1, kRed, kMagenta, kBlack};

   // --- Page 1: kinematic curves with chi2 cuts ---
   TCanvas *c1 = new TCanvas("c1", "Kinematics with chi2 cuts", 1800, 600);
   c1->Divide(3, 1);

   // Truth
   c1->cd(1);
   std::vector<double> tTh, tKE;
   for (auto &d : data) { tTh.push_back(d.thTrue); tKE.push_back(d.keTrue); }
   auto *gT = new TGraph(tTh.size(), tTh.data(), tKE.data());
   gT->SetTitle("MC Truth;#theta_{lab} [deg];KE [MeV]");
   gT->SetMarkerStyle(20); gT->SetMarkerSize(0.3); gT->SetMarkerColor(kBlack);
   gT->Draw("AP");

   // All reco
   c1->cd(2);
   std::vector<double> rTh, rKE;
   for (auto &d : data) { rTh.push_back(d.thReco); rKE.push_back(d.keReco); }
   auto *gR = new TGraph(rTh.size(), rTh.data(), rKE.data());
   gR->SetTitle(Form("All reco (N=%d);#theta_{lab} [deg];KE [MeV]", (int)data.size()));
   gR->SetMarkerStyle(20); gR->SetMarkerSize(0.3); gR->SetMarkerColor(kBlack);
   gR->Draw("AP");

   // Best cut overlay
   c1->cd(3);
   auto *gT2 = new TGraph(tTh.size(), tTh.data(), tKE.data());
   gT2->SetTitle("#chi^{2}/ndf < 0.2 (blue) vs all (gray);#theta_{lab} [deg];KE [MeV]");
   gT2->SetMarkerStyle(20); gT2->SetMarkerSize(0.2); gT2->SetMarkerColor(kGray);
   gT2->Draw("AP");
   std::vector<double> cTh, cKE;
   for (auto &d : data) { if (d.chi2 < 0.2) { cTh.push_back(d.thReco); cKE.push_back(d.keReco); } }
   auto *gC = new TGraph(cTh.size(), cTh.data(), cKE.data());
   gC->SetMarkerStyle(20); gC->SetMarkerSize(0.3); gC->SetMarkerColor(kBlue);
   gC->Draw("P SAME");
   auto *leg0 = new TLegend(0.5, 0.7, 0.88, 0.88);
   leg0->AddEntry(gT2, Form("Truth (%d)", (int)tTh.size()), "p");
   leg0->AddEntry(gC, Form("#chi^{2}<0.2 (%d)", (int)cTh.size()), "p");
   leg0->Draw();

   c1->SaveAs("data/chi2_cuts_kinematics.png");

   // --- Page 2: KE error distributions for each cut ---
   TCanvas *c2 = new TCanvas("c2", "KE error with chi2 cuts", 1800, 600);
   c2->Divide(3, 2);

   for (int ic = 0; ic < nCuts; ic++) {
      c2->cd(ic + 1);
      TH1F *h = new TH1F(Form("hKE_%d", ic),
                          Form("#chi^{2}/ndf < %.1f;(KE_{reco}-KE_{true})/KE_{true} [%%];Events", cuts[ic]),
                          80, -30, 30);
      int n = 0;
      for (auto &d : data) {
         if (d.chi2 < cuts[ic]) { h->Fill(d.keErr); n++; }
      }
      h->SetLineColor(colors[ic]);
      h->SetLineWidth(2);
      h->Fit("gaus", "Q");
      h->Draw("hist");
      if (h->GetFunction("gaus")) h->GetFunction("gaus")->Draw("same");

      double mean = h->GetMean();
      double rms = h->GetStdDev();
      TLatex *txt = new TLatex();
      txt->SetNDC(); txt->SetTextSize(0.04);
      txt->DrawLatex(0.15, 0.85, Form("N = %d (%.0f%%)", n, 100.0*n/data.size()));
      txt->DrawLatex(0.15, 0.80, Form("Mean = %.2f%%", mean));
      txt->DrawLatex(0.15, 0.75, Form("RMS = %.2f%%", rms));
   }

   c2->SaveAs("data/chi2_cuts_KE_error.png");

   // --- Page 3: summary table ---
   std::cout << "\n=== Chi2 cut summary ===" << std::endl;
   std::cout << "Cut       N      Eff%    KE_mean%  KE_RMS%  Th_mean   Th_RMS" << std::endl;
   for (int ic = 0; ic < nCuts; ic++) {
      int n = 0;
      double sumKE = 0, sumKE2 = 0, sumTh = 0, sumTh2 = 0;
      for (auto &d : data) {
         if (d.chi2 < cuts[ic]) {
            n++;
            double thErr = d.thReco - d.thTrue;
            sumKE += d.keErr; sumKE2 += d.keErr * d.keErr;
            sumTh += thErr; sumTh2 += thErr * thErr;
         }
      }
      if (n == 0) continue;
      double mKE = sumKE / n, rKE = std::sqrt(sumKE2 / n - mKE * mKE);
      double mTh = sumTh / n, rTh = std::sqrt(sumTh2 / n - mTh * mTh);
      printf("< %.1f   %4d   %5.1f%%   %+6.2f    %5.2f    %+5.2f     %.2f\n",
             cuts[ic], n, 100.0 * n / data.size(), mKE, rKE, mTh, rTh);
   }

   std::cout << "\nSaved: data/chi2_cuts_kinematics.png" << std::endl;
   std::cout << "Saved: data/chi2_cuts_KE_error.png" << std::endl;
}
