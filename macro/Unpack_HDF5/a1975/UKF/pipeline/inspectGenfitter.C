/// @file inspectGenfitter.C
/// @brief First-look diagnostics for an AtGenfitter output file.
/// Dumps per-event yields + fit quality, and writes diagnostic PNGs.
///   root -b -q 'pipeline/inspectGenfitter.C("/tmp/run_0106_genfitter_val.root", 30)'
/// nDump = number of per-track lines to print to screen (0 = none).

void inspectGenfitter(TString file = "/tmp/run_0106_genfitter_val.root", Int_t nDump = 30)
{
   gSystem->Load("libAtReconstruction.so");
   TFile f(file);
   TTree *t = (TTree *)f.Get("cbmsim");
   if (!t) { std::cout << "no cbmsim tree\n"; return; }

   TClonesArray *trkEvt = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &trkEvt);

   // histograms
   TH1F *hKE = new TH1F("hKE", "Fitted deuteron KE;KE [MeV];tracks", 100, 0, 60);
   TH1F *hTh = new TH1F("hTh", "Fitted #theta;#theta [deg];tracks", 90, 0, 180);
   TH1F *hChi2ndf = new TH1F("hChi2ndf", "#chi^{2}/ndf;#chi^{2}/ndf;tracks", 100, 0, 50);
   TH1F *hVz = new TH1F("hVz", "Vertex z;z [mm];tracks", 120, -100, 1100);
   TH1F *hVr = new TH1F("hVr", "Vertex r=#sqrt{x^{2}+y^{2}};r [mm];tracks", 100, 0, 80);
   TH2F *hKEth = new TH2F("hKEth", "KE vs #theta;#theta [deg];KE [MeV]", 90, 0, 180, 100, 0, 60);

   Long64_t N = t->GetEntries();
   long nEvtWithTrk = 0, nFitted = 0, nConv = 0, nPrinted = 0;
   // population counters
   long nKEzero = 0, nKEhuge = 0, nPerp = 0, nFwd = 0, nFwdConv = 0;
   double fwdKEsum = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      auto *ev = (AtTrackingEvent *)trkEvt->At(0);
      if (!ev) continue;
      auto &fitted = ev->GetFittedTracks();
      if (!fitted.empty()) nEvtWithTrk++;
      for (auto &ftp : fitted) {
         auto *ft = ftp.get();
         if (!ft) continue;
         nFitted++;
         const auto &kin = ft->GetKinematics(0);
         const auto &vtx = ft->GetVertex(0);
         auto &meta = ft->GetTrackMetadata();
         double chi2 = meta ? meta->GetChi2() : -1;
         int ndf = meta ? meta->GetNdf() : 0;
         bool conv = meta ? meta->GetFitConverged() : false;
         double c2n = (ndf > 0) ? chi2 / ndf : -1;
         if (conv) nConv++;

         double thdeg = kin.theta * TMath::RadToDeg();
         double ke = kin.kineticEnergy;
         if (ke < 0.5) nKEzero++;
         if (ke > 60) nKEhuge++;
         if (thdeg > 70 && thdeg < 110) nPerp++;
         if (thdeg < 45) { nFwd++; if (conv && ke > 0.5 && ke < 60) { nFwdConv++; fwdKEsum += ke; } }
         hKE->Fill(kin.kineticEnergy);
         hTh->Fill(thdeg);
         if (c2n >= 0) hChi2ndf->Fill(c2n);
         hVz->Fill(vtx.Z());
         hVr->Fill(std::hypot(vtx.X(), vtx.Y()));
         hKEth->Fill(thdeg, kin.kineticEnergy);

         if (nPrinted < nDump) {
            printf("evt %5lld trk %2d  KE=%7.2f MeV  th=%6.2f deg  phi=%7.2f  vtx=(%6.1f,%6.1f,%6.1f)mm  "
                   "chi2/ndf=%7.2f/%d=%6.2f  conv=%d\n",
                   i, ft->GetTrackID(), kin.kineticEnergy, thdeg, kin.phi * TMath::RadToDeg(), vtx.X(), vtx.Y(),
                   vtx.Z(), chi2, ndf, c2n, conv);
            nPrinted++;
         }
      }
   }

   printf("\n=== SUMMARY (%lld events) ===\n", N);
   printf("events with >=1 fitted track : %ld\n", nEvtWithTrk);
   printf("total fitted tracks          : %ld\n", nFitted);
   printf("converged (genfit flag)      : %ld (%.1f%%)\n", nConv, 100.0 * nConv / std::max(1L, nFitted));
   printf("KE     mean=%.2f rms=%.2f MeV\n", hKE->GetMean(), hKE->GetRMS());
   printf("theta  mean=%.2f rms=%.2f deg\n", hTh->GetMean(), hTh->GetRMS());
   printf("chi2/ndf mean=%.2f median~=%.2f\n", hChi2ndf->GetMean(),
          hChi2ndf->GetBinCenter(hChi2ndf->GetMaximumBin()));
   printf("vtx z  mean=%.1f rms=%.1f mm | vtx r mean=%.1f rms=%.1f mm\n", hVz->GetMean(), hVz->GetRMS(),
          hVr->GetMean(), hVr->GetRMS());
   printf("\n--- populations (of %ld fitted) ---\n", nFitted);
   printf("KE < 0.5 MeV  (collapsed) : %ld (%.1f%%)\n", nKEzero, 100.0 * nKEzero / nFitted);
   printf("KE > 60  MeV  (diverged)  : %ld (%.1f%%)\n", nKEhuge, 100.0 * nKEhuge / nFitted);
   printf("theta 70-110 (perp blob)  : %ld (%.1f%%)\n", nPerp, 100.0 * nPerp / nFitted);
   printf("theta < 45   (forward)    : %ld (%.1f%%)\n", nFwd, 100.0 * nFwd / nFitted);
   printf("  forward & conv & physical KE: %ld  <KE>=%.2f MeV\n", nFwdConv, nFwdConv ? fwdKEsum / nFwdConv : 0);

   gStyle->SetOptStat(1110);
   TCanvas *c = new TCanvas("c", "genfitter", 1500, 1000);
   c->Divide(3, 2);
   c->cd(1); hKE->Draw();
   c->cd(2); hTh->Draw();
   c->cd(3); hChi2ndf->Draw();
   c->cd(4); hVz->Draw();
   c->cd(5); hVr->Draw();
   c->cd(6); hKEth->Draw("colz");
   TString png = file; png.ReplaceAll(".root", "_inspect.png");
   c->SaveAs(png);
   printf("\nwrote %s\n", png.Data());
}
