/// @file phi_bias_diagnostic.C
/// @brief Compare UKF φ to MC φ at TWO reference points:
///   (i)  φ_MC at the production vertex (what we did before)
///   (ii) φ_MC at the FIRST AtMCPoint inside drift_volume (the spatial
///        location where the UKF actually anchors its angle)
///
/// If (i) - UKF shows a positive bias and (ii) - UKF is ~0, the bias is
/// just the helix's φ rotation between vertex and first cluster.

void phi_bias_diagnostic(int P = 800)
{
   const double mass_pi = 139.57039;
   const double radToMrad = 1000.;

   TFile fS(Form("data/HYDRAsim_p%d.root", P));
   TFile fU(Form("data/output_ukf_HYDRA_p%d.root", P));
   auto *tS = (TTree *)fS.Get("cbmsim");
   auto *tU = (TTree *)fU.Get("cbmsim");
   auto *trks = new TClonesArray("AtMCTrack");
   auto *pts  = new TClonesArray("AtMCPoint");
   auto *te   = new TClonesArray("AtTrackingEvent");
   tS->SetBranchAddress("MCTrack", &trks);
   tS->SetBranchAddress("AtTpcPoint", &pts);
   tU->SetBranchAddress("AtTrackingEvent", &te);

   TH1F h1("h1", Form("p=%d MeV/c;#Delta#varphi (mrad);entries", P), 80, -300, 300);
   TH1F h2("h2", Form("p=%d MeV/c;#Delta#varphi (mrad);entries", P), 80, -300, 300);
   h1.SetDirectory(nullptr); h2.SetDirectory(nullptr);
   h1.SetLineColor(kRed + 1);   h1.SetLineWidth(2);
   h2.SetLineColor(kBlue + 1);  h2.SetLineWidth(2);

   Long64_t n = std::min(tS->GetEntries(), tU->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tS->GetEntry(i);
      tU->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double phMCvertex = std::atan2(mc->GetPy(), mc->GetPx());

      // Find primary-pion track id and the MC point with smallest x in drift_volume
      double xMin = 1e9;
      double phMCfirst = -999;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (p->GetTrackID() != 0) continue;
         double x_mm = p->GetX() * 10.;
         if (x_mm < xMin) {
            xMin = x_mm;
            phMCfirst = std::atan2(p->GetPy(), p->GetPx());
         }
      }
      if (phMCfirst == -999) continue;

      if (te->GetEntries() == 0) continue;
      auto *trkEvt = (AtTrackingEvent *)te->At(0);
      auto &fitted = trkEvt->GetFittedTracks();
      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double c = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
         if (c < bestChi) { bestChi = c; best = t.get(); }
      }
      if (!best) continue;
      double phUKF = best->GetKinematics().phi;
      auto wrap = [](double d) {
         while (d > M_PI) d -= 2 * M_PI;
         while (d < -M_PI) d += 2 * M_PI;
         return d;
      };
      h1.Fill(wrap(phUKF - phMCvertex) * radToMrad);
      h2.Fill(wrap(phUKF - phMCfirst)  * radToMrad);
   }

   auto fitCore = [](TH1F &h, double &mu, double &sig) {
      double rms = h.GetRMS();
      TFitResultPtr fr = h.Fit("gaus", "SQR0", "", -3 * rms, 3 * rms);
      mu = fr.Get() ? fr->Parameter(1) : h.GetMean();
      sig = fr.Get() ? fr->Parameter(2) : rms;
   };
   double mu1, sig1, mu2, sig2;
   fitCore(h1, mu1, sig1);
   fitCore(h2, mu2, sig2);

   std::cout << "\n=== p=" << P << " MeV/c ===\n";
   std::cout << "  (i)  UKF - MC@vertex     : bias = " << mu1 << " mrad, σ = " << sig1 << " mrad\n";
   std::cout << "  (ii) UKF - MC@firstMCpt  : bias = " << mu2 << " mrad, σ = " << sig2 << " mrad\n";

   auto *c = new TCanvas(Form("c_%d", P), "", 800, 600);
   c->SetGrid();
   h1.Draw("HIST");
   h2.Draw("HIST same");
   h1.SetMaximum(1.3 * std::max(h1.GetMaximum(), h2.GetMaximum()));
   auto *leg = new TLegend(0.6, 0.75, 0.93, 0.90);
   leg->SetBorderSize(0); leg->SetFillStyle(0);
   leg->AddEntry(&h1, Form("UKF #minus MC@vertex (#mu=%.1f mrad)", mu1), "l");
   leg->AddEntry(&h2, Form("UKF #minus MC@firstMC (#mu=%.1f mrad)", mu2), "l");
   leg->Draw();
   c->SaveAs(Form("data/phi_bias_p%d.png", P));
}
