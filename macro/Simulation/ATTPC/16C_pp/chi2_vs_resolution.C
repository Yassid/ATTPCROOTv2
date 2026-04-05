void chi2_vs_resolution()
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
   std::vector<double> vChi2, vKEerr, vAbsKEerr, vThErr, vAbsThErr, vVtxR;

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
      double keR = kin.kineticEnergy;
      double thR = kin.theta * 180.0 / M_PI;
      if (keR < 0.05 || keR > 50 || thR < 10 || thR > 170) continue;

      auto [pval, chi2, bchi2, ndf, bndf, conv] = ft[0]->GetStats();
      if (ndf <= 0 || chi2 <= 0) continue;

      // Initial vertex distance
      if (!pa || pa->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)pa->At(0);
      auto &tracks = pe->GetTrackCand();
      if (tracks.empty()) continue;
      int best = 0;
      for (size_t t = 1; t < tracks.size(); t++)
         if (tracks[t].GetHitArray().size() > tracks[best].GetHitArray().size()) best = t;
      auto *cl = tracks[best].GetHitClusterArray();
      if (cl->empty()) continue;
      auto pos0 = cl->front().GetPosition();
      double vtxR = std::sqrt(pos0.X() * pos0.X() + pos0.Y() * pos0.Y());

      double keErr = (keR - keT) / keT * 100;
      double thErr = thR - thT;

      vChi2.push_back(chi2);
      vKEerr.push_back(keErr);
      vAbsKEerr.push_back(std::abs(keErr));
      vThErr.push_back(thErr);
      vAbsThErr.push_back(std::abs(thErr));
      vVtxR.push_back(vtxR);
   }

   std::cout << "Events: " << vChi2.size() << std::endl;

   TCanvas *c = new TCanvas("c", "chi2 vs resolution", 1400, 900);
   c->Divide(3, 2);

   // chi2 vs KE error
   c->cd(1);
   auto *h1 = new TH2F("h1", "#chi^{2}/ndf vs KE error;#chi^{2}/ndf;(KE_{reco}-KE_{true})/KE_{true} [%]",
                        50, 0, 1.0, 100, -30, 30);
   for (size_t i = 0; i < vChi2.size(); i++) h1->Fill(vChi2[i], vKEerr[i]);
   h1->Draw("COLZ");

   // chi2 vs theta error
   c->cd(2);
   auto *h2 = new TH2F("h2", "#chi^{2}/ndf vs #theta error;#chi^{2}/ndf;#theta_{reco}-#theta_{true} [deg]",
                        50, 0, 1.0, 100, -5, 5);
   for (size_t i = 0; i < vChi2.size(); i++) h2->Fill(vChi2[i], vThErr[i]);
   h2->Draw("COLZ");

   // chi2 vs vtxR
   c->cd(3);
   auto *h3 = new TH2F("h3", "#chi^{2}/ndf vs vtxR;#chi^{2}/ndf;vtxR [mm]",
                        50, 0, 1.0, 50, 0, 50);
   for (size_t i = 0; i < vChi2.size(); i++) h3->Fill(vChi2[i], vVtxR[i]);
   h3->Draw("COLZ");

   // |KE error| profile vs chi2
   c->cd(4);
   auto *h4 = new TH2F("h4", "|KE error| vs #chi^{2}/ndf;#chi^{2}/ndf;|KE error| [%]",
                        25, 0, 1.0, 50, 0, 20);
   for (size_t i = 0; i < vChi2.size(); i++) h4->Fill(vChi2[i], vAbsKEerr[i]);
   auto *p4 = h4->ProfileX("p4");
   p4->SetTitle("mean |KE error| vs #chi^{2}/ndf;#chi^{2}/ndf;mean |KE error| [%]");
   p4->SetMarkerStyle(20);
   p4->SetMarkerSize(0.6);
   p4->SetLineColor(kBlue);
   p4->Draw("E1");

   // |theta error| profile vs chi2
   c->cd(5);
   auto *h5 = new TH2F("h5", "|#theta error| vs #chi^{2}/ndf;#chi^{2}/ndf;|#theta error| [deg]",
                        25, 0, 1.0, 50, 0, 5);
   for (size_t i = 0; i < vChi2.size(); i++) h5->Fill(vChi2[i], vAbsThErr[i]);
   auto *p5 = h5->ProfileX("p5");
   p5->SetTitle("mean |#theta error| vs #chi^{2}/ndf;#chi^{2}/ndf;mean |#theta error| [deg]");
   p5->SetMarkerStyle(20);
   p5->SetMarkerSize(0.6);
   p5->SetLineColor(kRed);
   p5->Draw("E1");

   // KE resolution in chi2 bins: RMS of KE error for low vs high chi2
   c->cd(6);
   TH1F *hLo = new TH1F("hLo", "#chi^{2}/ndf < 0.1", 60, -30, 30);
   TH1F *hHi = new TH1F("hHi", "#chi^{2}/ndf > 0.3", 60, -30, 30);
   for (size_t i = 0; i < vChi2.size(); i++) {
      if (vChi2[i] < 0.1) hLo->Fill(vKEerr[i]);
      else if (vChi2[i] > 0.3) hHi->Fill(vKEerr[i]);
   }
   hLo->SetLineColor(kBlue);
   hHi->SetLineColor(kRed);
   hLo->SetTitle(Form("KE error: low #chi^{2} (RMS=%.1f%%, N=%d) vs high (RMS=%.1f%%, N=%d);KE error [%%];Events",
                       hLo->GetStdDev(), (int)hLo->GetEntries(), hHi->GetStdDev(), (int)hHi->GetEntries()));
   hLo->Draw("hist");
   hHi->Draw("hist same");
   auto *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg->AddEntry(hLo, Form("#chi^{2}<0.1 (N=%d)", (int)hLo->GetEntries()), "l");
   leg->AddEntry(hHi, Form("#chi^{2}>0.3 (N=%d)", (int)hHi->GetEntries()), "l");
   leg->Draw();

   c->SaveAs("data/chi2_vs_resolution.png");
   std::cout << "Saved: data/chi2_vs_resolution.png" << std::endl;
}
