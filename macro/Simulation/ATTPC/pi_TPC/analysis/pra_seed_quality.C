/// @brief Quantify PRA seed-parameter quality (R, theta, phi) vs MC truth.
/// Useful to identify which PRA seed component limits UKF resolution.
void pra_seed_quality()
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tDigi->SetBranchAddress("AtPatternEvent", &patArr);

   const double Bz_T = 0.5;
   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<TH1F *> hdR(NB), hdTh(NB), hdPhi(NB);
   for (int b = 0; b < NB; ++b) {
      hdR[b] = new TH1F(Form("hdR_%d", b), Form("th %.0f-%.0f;dR/R_truth (%%);counts", edges[b], edges[b + 1]),
                        80, -50., 50.);
      hdTh[b] = new TH1F(Form("hdTh_%d", b), Form(";dtheta (deg);counts"), 80, -30., 30.);
      hdPhi[b] = new TH1F(Form("hdPhi_%d", b), Form(";dphi (deg);counts"), 80, -30., 30.);
   }

   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double pmc = std::sqrt(mc->GetPx()*mc->GetPx() + mc->GetPy()*mc->GetPy() + mc->GetPz()*mc->GetPz());
      double pmc_MeV = pmc * 1000.;
      double thMC_lab = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      double phMC = std::atan2(mc->GetPy(), mc->GetPx()) * 180. / M_PI;
      double pT_MC = pmc_MeV * std::sin(thMC_lab * M_PI / 180.);
      double R_MC_mm = pT_MC / (0.3 * Bz_T);

      int b = -1;
      for (int k = 0; k < NB; ++k) if (thMC_lab >= edges[k] && thMC_lab < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;
      if (patArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patArr->At(0);
      auto &cands = pat->GetTrackCand();
      if (cands.empty()) continue;
      auto &t = cands[0];
      double R_PRA = t.GetGeoRadius();
      double thPRA_digi = t.GetGeoTheta() * 180. / M_PI;
      double phPRA = t.GetGeoPhi() * 180. / M_PI;
      double thPRA_lab = 180. - thPRA_digi;
      if (!std::isfinite(R_PRA) || R_PRA <= 0 || !std::isfinite(thPRA_digi)) continue;

      double dR_pct = (R_PRA - R_MC_mm) / R_MC_mm * 100.;
      double dTh = thPRA_lab - thMC_lab;
      double dPhi = phPRA - phMC;
      while (dPhi > 180) dPhi -= 360;
      while (dPhi < -180) dPhi += 360;

      hdR[b]->Fill(dR_pct);
      hdTh[b]->Fill(dTh);
      hdPhi[b]->Fill(dPhi);
   }

   std::cout << "\n=== PRA seed quality vs MC truth (pi+ isotropic) ===\n";
   std::cout << "theta_MC      N    dR (Gauss)         dTheta_lab (Gauss)    dPhi (Gauss)\n";
   std::cout << std::string(78, '-') << "\n";
   auto fitG = [](TH1F *h, double lo, double hi) -> std::tuple<double, double, double> {
      if (h->GetEntries() < 15) return {NAN, NAN, h->GetEntries()};
      double rms = h->GetRMS();
      double mean = h->GetMean();
      double a = std::max(lo, mean - 2.5 * rms);
      double bb = std::min(hi, mean + 2.5 * rms);
      h->Fit("gaus", "Q0", "", a, bb);
      auto *f = h->GetFunction("gaus");
      if (!f) return {mean, rms, h->GetEntries()};
      return {f->GetParameter(1), f->GetParameter(2), h->GetEntries()};
   };
   for (int b = 0; b < NB; ++b) {
      auto [muR, sgR, nR] = fitG(hdR[b], -50, 50);
      auto [muT, sgT, nT] = fitG(hdTh[b], -30, 30);
      auto [muP, sgP, nP] = fitG(hdPhi[b], -30, 30);
      char buf[256];
      snprintf(buf, sizeof(buf),
               "%4.0f-%4.0f  %4d   %+5.1f±%4.1f%%      %+5.2f±%4.2f°       %+5.2f±%4.2f°",
               edges[b], edges[b + 1], (int)nR, muR, sgR, muT, sgT, muP, sgP);
      std::cout << buf << "\n";
   }
}
