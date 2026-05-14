/// @file quick_perf.C
/// @brief Quick σ_KE + fit-yield summary from output_reco_ukf.root vs attpcsim.root

void quick_perf(const char *fitF = "data/output_reco_ukf.root",
                const char *simF = "data/attpcsim.root")
{
   TFile fS(simF), fF(fitF);
   auto *tS = (TTree *)fS.Get("cbmsim");
   auto *tF = (TTree *)fF.Get("cbmsim");
   if (!tS || !tF) { std::cerr << "missing trees\n"; return; }

   auto *trks = new TClonesArray("AtMCTrack");
   auto *te   = new TClonesArray("AtTrackingEvent");
   tS->SetBranchAddress("MCTrack", &trks);
   tF->SetBranchAddress("AtTrackingEvent", &te);

   const double mp = 938.272;
   Long64_t n = std::min(tS->GetEntries(), tF->GetEntries());
   int nFit = 0, nGood = 0;
   double sumChi = 0, sumDKE = 0, sumDKE2 = 0;

   auto *h = new TH1F("hDKE", ";(KE_{fit}-KE_{MC})/KE_{MC};", 100, -0.5, 0.5);

   for (Long64_t i = 0; i < n; ++i) {
      tS->GetEntry(i); tF->GetEntry(i);
      if (trks->GetEntries() < 2) continue;
      // proton is the second MC track in 16C+p
      AtMCTrack *mc = nullptr;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *t = (AtMCTrack *)trks->At(j);
         if (t->GetPdgCode() == 2212) { mc = t; break; }
      }
      if (!mc) continue;
      double Px = mc->GetPx(), Py = mc->GetPy(), Pz = mc->GetPz();
      double pmc = std::sqrt(Px*Px + Py*Py + Pz*Pz) * 1000.;
      double keMC = std::sqrt(pmc*pmc + mp*mp) - mp;

      if (te->GetEntries() == 0) continue;
      auto *trkEvt = (AtTrackingEvent *)te->At(0);
      auto &fitted = trkEvt->GetFittedTracks();
      AtFittedTrack *best = nullptr; double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2()/ndf : 1e30;
         if (cc < bestChi) { bestChi = cc; best = t.get(); }
      }
      if (!best) continue;
      nFit++;
      double keFit = best->GetKinematics().kineticEnergy;
      double keXtr = best->GetKinematicsXtr().kineticEnergy;
      double dke = (keFit - keMC) / keMC;
      if (nFit <= 20) printf("evt %lld: keMC=%.2f  keXtr=%.2f  keFit(@vtx)=%.2f  chi2/ndf=%.2f\n",
                              i, keMC, keXtr, keFit, bestChi);
      sumChi += bestChi;
      sumDKE += dke;
      sumDKE2 += dke*dke;
      h->Fill(dke);
      if (std::abs(dke) < 0.1) nGood++;
   }

   double mean = nFit ? sumDKE/nFit : 0;
   double rms  = nFit ? std::sqrt(std::max(0., sumDKE2/nFit - mean*mean)) : 0;
   printf("Events checked: %lld\n", n);
   printf("Fitted:         %d  (%.1f%%)\n", nFit, 100.*nFit/n);
   printf("|ΔKE/KE|<0.1:   %d  (%.1f%% of fit)\n", nGood, nFit?100.*nGood/nFit:0);
   printf("ΔKE/KE  mean = %+.3f  RMS = %.3f\n", mean, rms);
   printf("⟨chi2/ndf⟩    = %.3f\n", nFit ? sumChi/nFit : 0);

   h->Fit("gaus", "Q0");
   auto *g = h->GetFunction("gaus");
   if (g) printf("Gaussian core:  μ = %+.3f  σ = %.3f\n", g->GetParameter(1), g->GetParameter(2));
}
