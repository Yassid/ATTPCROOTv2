void vtxR_vs_error2()
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
   std::vector<double> vR, vKEerr, vThErr, vKEtrue;

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

      if (!fa || fa->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)fa->At(0);
      auto &ft = te->GetFittedTracks();
      if (ft.empty()) continue;
      auto kin = ft[0]->GetKinematics();
      double keR = kin.kineticEnergy;
      double thR = kin.theta * 180.0 / M_PI;
      if (keR < 0.05 || keR > 50 || thR < 10 || thR > 170) continue;

      vR.push_back(vtxR);
      vKEerr.push_back((keR - keT) / keT * 100);
      vThErr.push_back(thR - thT);
      vKEtrue.push_back(keT);
   }

   TCanvas *c = new TCanvas("c", "vtxR vs errors", 1400, 900);
   c->Divide(2, 2);

   // KE error vs vtxR (zoomed)
   c->cd(1);
   auto *h1 = new TH2F("f1", "KE error vs initial vtxR;vtxR (first cluster) [mm];(KE_{reco}-KE_{true})/KE_{true} [%]",
                        100, 0, 50, 100, -30, 30);
   for (size_t i = 0; i < vR.size(); i++) h1->Fill(vR[i], vKEerr[i]);
   h1->Draw("COLZ");

   // Theta error vs vtxR (zoomed)
   c->cd(2);
   auto *h2 = new TH2F("f2", "#theta error vs initial vtxR;vtxR (first cluster) [mm];#theta_{reco}-#theta_{true} [deg]",
                        100, 0, 50, 100, -5, 5);
   for (size_t i = 0; i < vR.size(); i++) h2->Fill(vR[i], vThErr[i]);
   h2->Draw("COLZ");

   // |KE error| profile
   c->cd(3);
   auto *h3 = new TH2F("f3", "|KE error| vs vtxR;vtxR [mm];|KE error| [%]", 25, 0, 50, 50, 0, 20);
   for (size_t i = 0; i < vR.size(); i++) h3->Fill(vR[i], std::abs(vKEerr[i]));
   auto *p3 = h3->ProfileX("p3");
   p3->SetTitle("|KE error| profile vs initial vtxR;vtxR (first cluster) [mm];mean |KE error| [%]");
   p3->SetMarkerStyle(20);
   p3->SetMarkerSize(0.6);
   p3->SetLineColor(kBlue);
   p3->Draw("E1");

   // vtxR vs KE_true
   c->cd(4);
   auto *h4 = new TH2F("f4", "vtxR vs KE_{true};KE_{true} [MeV];vtxR (first cluster) [mm]",
                        100, 0, 20, 100, 0, 50);
   for (size_t i = 0; i < vR.size(); i++) h4->Fill(vKEtrue[i], vR[i]);
   h4->Draw("COLZ");

   c->SaveAs("data/vtxR_vs_error.png");
   std::cout << "Saved: data/vtxR_vs_error.png" << std::endl;
}
