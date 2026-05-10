/// @brief KE resolution and theta resolution as a function of theta_MC.
///
/// Reads ./data/attpcsim.root + ./data/output_ukf_only.root, bins fits by
/// theta_MC, fits a Gaussian core to (KE_fit - KE_MC) per bin, and tabulates
/// fit yield, KE bias, sigma, sigma/<KE>, and theta-residual sigma.
///
/// Run: root -b -q ke_vs_theta.C

void ke_vs_theta(const char *ukfFile = "data/output_ukf_only.root", const char *tag = "PRA")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TFile fSim("data/attpcsim.root");
   TFile fUKF(ukfFile);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<TH1F *> hKE(NB), hTh(NB);
   std::vector<int> nThrown(NB, 0);
   std::vector<double> sumKE(NB, 0.0);

   for (int b = 0; b < NB; ++b) {
      hKE[b] = new TH1F(Form("hKE_%d", b),
                        Form("theta %.0f-%.0f deg;dKE (MeV);counts", edges[b], edges[b + 1]),
                        80, -25., 25.);
      hTh[b] = new TH1F(Form("hTh_%d", b),
                        Form("theta %.0f-%.0f deg;dtheta (deg);counts", edges[b], edges[b + 1]),
                        80, -10., 10.);
   }

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   int nThrownTot = 0, nFitTot = 0;
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;

      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
      double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      ++nThrownTot;
      int bin = -1;
      for (int b = 0; b < NB; ++b) {
         if (thMC >= edges[b] && thMC < edges[b + 1]) { bin = b; break; }
      }
      if (bin < 0) continue;
      ++nThrown[bin];
      sumKE[bin] += KEmc;

      if (teArr->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      auto &fitted = te->GetFittedTracks();
      if (fitted.empty()) continue;

      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double chi = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
         if (chi < bestChi) { bestChi = chi; best = t.get(); }
      }
      if (!best) continue;
      double KEfit = best->GetKinematics().kineticEnergy;
      double thFit = best->GetKinematics().theta * 180. / M_PI;
      hKE[bin]->Fill(KEfit - KEmc);
      hTh[bin]->Fill(thFit - thMC);
      ++nFitTot;
   }

   std::cout << "\n=== KE / theta resolution vs theta_MC (isotropic pi+, B=0.5T, P10@1bar, 5-50 MeV) ===\n";
   std::cout << "Thrown: " << nThrownTot << "    Fitted: " << nFitTot
             << "    Yield: " << 100. * nFitTot / std::max(1, nThrownTot) << " %\n\n";
   std::cout << std::left
             << std::setw(15) << "theta_MC (deg)"
             << std::setw(8)  << "Nthrown"
             << std::setw(7)  << "Nfit"
             << std::setw(7)  << "yield"
             << std::setw(11) << "<KE>(MeV)"
             << std::setw(13) << "KE_bias"
             << std::setw(13) << "KE_sigma"
             << std::setw(11) << "sigma/E"
             << std::setw(13) << "th_sigma"
             << "\n";
   std::cout << std::string(98, '-') << "\n";

   auto fitGaus = [](TH1F *h, double lo, double hi) -> std::pair<double, double> {
      if (h->GetEntries() < 15) return {std::nan(""), std::nan("")};
      double rms = h->GetRMS();
      double mean = h->GetMean();
      double a = std::max(lo, mean - 2.5 * rms);
      double b = std::min(hi, mean + 2.5 * rms);
      h->Fit("gaus", "Q0", "", a, b);
      auto *f = h->GetFunction("gaus");
      if (!f) return {mean, rms};
      return {f->GetParameter(1), f->GetParameter(2)};
   };

   for (int b = 0; b < NB; ++b) {
      auto [muKE, sgKE] = fitGaus(hKE[b], -25., 25.);
      auto [muTh, sgTh] = fitGaus(hTh[b], -10., 10.);
      double meanKE = nThrown[b] > 0 ? sumKE[b] / nThrown[b] : 0;
      double yld = nThrown[b] > 0 ? 100. * hKE[b]->GetEntries() / nThrown[b] : 0;
      double rel = (meanKE > 0 && std::isfinite(sgKE)) ? sgKE / meanKE * 100. : 0;
      char buf[256];
      snprintf(buf, sizeof(buf),
               "%-15s %-8d %-7d %5.1f%%  %8.2f    bias=%+5.2f   sg=%5.2f    %5.1f%%   sg=%5.2f",
               Form("%4.0f-%4.0f", edges[b], edges[b + 1]),
               nThrown[b], (int)hKE[b]->GetEntries(), yld, meanKE,
               muKE, sgKE, rel, sgTh);
      std::cout << buf << "\n";
   }

   // Plot sigma/E vs theta
   auto *c = new TCanvas("c", "KE/theta resolution vs theta_MC", 1400, 600);
   c->Divide(2, 1);
   const int NP = NB;
   double thMid[NP], thErr[NP], sgRel[NP], sgRelErr[NP], sgThArr[NP], sgThErr[NP];
   for (int b = 0; b < NB; ++b) {
      thMid[b] = 0.5 * (edges[b] + edges[b + 1]);
      thErr[b] = 0.5 * (edges[b + 1] - edges[b]);
      auto [muKE, sgKE] = fitGaus(hKE[b], -25., 25.);
      auto [muTh, sgTh] = fitGaus(hTh[b], -10., 10.);
      double meanKE = nThrown[b] > 0 ? sumKE[b] / nThrown[b] : 1;
      sgRel[b] = (std::isfinite(sgKE) && meanKE > 0) ? sgKE / meanKE * 100. : 0;
      sgRelErr[b] = sgRel[b] / std::sqrt(2. * std::max(1., (double)hKE[b]->GetEntries()));
      sgThArr[b] = std::isfinite(sgTh) ? sgTh : 0;
      sgThErr[b] = sgThArr[b] / std::sqrt(2. * std::max(1., (double)hTh[b]->GetEntries()));
   }
   c->cd(1);
   gPad->SetGrid();
   auto *gKE = new TGraphErrors(NP, thMid, sgRel, thErr, sgRelErr);
   gKE->SetTitle("KE resolution vs theta_{MC};#theta_{MC} (deg);#sigma_{KE}/<KE> (%)");
   gKE->SetMarkerStyle(20);
   gKE->SetMarkerColor(kAzure + 2);
   gKE->SetLineColor(kAzure + 2);
   gKE->Draw("AP");
   c->cd(2);
   gPad->SetGrid();
   auto *gTh = new TGraphErrors(NP, thMid, sgThArr, thErr, sgThErr);
   gTh->SetTitle("theta resolution vs theta_{MC};#theta_{MC} (deg);#sigma_{#theta} (deg)");
   gTh->SetMarkerStyle(20);
   gTh->SetMarkerColor(kGreen + 2);
   gTh->SetLineColor(kGreen + 2);
   gTh->Draw("AP");
   c->SaveAs(Form("data/ke_vs_theta_%s.png", tag));
   c->SaveAs(Form("data/ke_vs_theta_%s.pdf", tag));
   std::cout << "\nWrote data/ke_vs_theta_" << tag << ".{png,pdf}\n";
}
