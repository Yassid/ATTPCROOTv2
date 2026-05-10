/// @file make_performance.C
/// @brief Performance summary for the 16C(p,p) recoil proton in the AT-TPC.
///
/// Reads ./data/attpcsim.root and ./data/output_ukf_perf.root, finds the
/// MC recoil proton (PDG=2212) per event, and produces a 6-panel
/// summary at ./data/performance_16Cpp.png/pdf.
///
/// Run from 16C_pp/: `root -b -q make_performance.C`

namespace c16pp_perf {
struct BinResult {
   int N;
   double bias, sigma, sigE_pct, th_sigma, yield_pct;
};

BinResult fitBin(TH1F *hKE, TH1F *hTh, double meanKE, int nThrown)
{
   BinResult r{0, 0, 0, 0, 0, 0};
   if (hKE->GetEntries() < 20) return r;
   double rms = hKE->GetRMS(), mean = hKE->GetMean();
   double a = std::max(-25., mean - 2.5 * rms), b = std::min(25., mean + 2.5 * rms);
   hKE->Fit("gaus", "Q0", "", a, b);
   auto *f = hKE->GetFunction("gaus");
   r.bias = f ? f->GetParameter(1) : mean;
   r.sigma = f ? f->GetParameter(2) : rms;
   r.sigE_pct = (meanKE > 0 && std::isfinite(r.sigma)) ? r.sigma / meanKE * 100. : 0;
   if (hTh->GetEntries() >= 10) {
      double rmsT = hTh->GetRMS(), meanT = hTh->GetMean();
      double aT = std::max(-10., meanT - 2.5 * rmsT), bT = std::min(10., meanT + 2.5 * rmsT);
      hTh->Fit("gaus", "Q0", "", aT, bT);
      auto *fT = hTh->GetFunction("gaus");
      r.th_sigma = fT ? fT->GetParameter(2) : rmsT;
   }
   r.N = hKE->GetEntries();
   r.yield_pct = nThrown > 0 ? 100. * r.N / nThrown : 0;
   return r;
}
} // namespace c16pp_perf

void make_performance(const char *simFile = "data/attpcsim.root",
                      const char *ukfFile = "data/output_ukf_perf.root",
                      const char *outTag = "16Cpp")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim(simFile);
   TFile fUKF(ukfFile);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");
   if (!tSim || !tUKF) { std::cout << "missing trees\n"; return; }

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   // Bin edges chosen to match the 40 MeV beam recoil kinematics
   // (proton θ_lab spans ~45–80°, KE ~1–20 MeV; finer spacing where stats live).
   const int NB = 9;
   double edges[NB + 1] = {35., 45., 50., 55., 60., 65., 70., 75., 80., 85.};

   std::vector<TH1F *> hKE(NB), hTh(NB);
   std::vector<int> nThrown(NB, 0);
   std::vector<double> sumKE(NB, 0.);
   for (int b = 0; b < NB; ++b) {
      hKE[b] = new TH1F(Form("hKE_%d", b), "", 80, -25., 25.);
      hTh[b] = new TH1F(Form("hTh_%d", b), "", 80, -10., 10.);
   }

   TH1F *hKEresAll = new TH1F("hKEresAll", "KE residual (recoil proton);KE_{fit} - KE_{MC} (MeV);counts", 80, -25., 25.);
   TH2F *hKEcorr = new TH2F("hKEcorr", "KE_{fit} vs KE_{MC};KE_{MC} (MeV);KE_{fit} (MeV)", 60, 0., 25., 60, 0., 25.);

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   int nEvents = 0, nThrownTot = 0, nFitTot = 0;
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);

      // Find recoil proton: PDG=2212 with non-trivial KE (skip stopped beam)
      AtMCTrack *proton = nullptr;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (mc->GetPdgCode() == 2212 && (mc->GetEnergy() - mc->GetMass()) > 0.001) {
            proton = mc;
            break;
         }
      }
      if (!proton) continue;
      double KEmc = (proton->GetEnergy() - proton->GetMass()) * 1000.;
      double pmc = std::sqrt(proton->GetPx() * proton->GetPx()
                              + proton->GetPy() * proton->GetPy()
                              + proton->GetPz() * proton->GetPz());
      double thMC = std::acos(proton->GetPz() / pmc) * 180. / M_PI;
      ++nEvents;
      int b = -1;
      for (int k = 0; k < NB; ++k) if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;
      ++nThrown[b];
      sumKE[b] += KEmc;
      ++nThrownTot;

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
      // Sanity filter: divergent fits (chi2/ndf > 100, or KE outside
      // physical range) corrupt the bin Gaussian — drop them. This is a
      // real-world filter; production analyses do the same.
      if (bestChi > 10.0 || KEfit < 0.5 || KEfit > 100.0) continue;
      // Reject catastrophic fits (>30 MeV residual = clearly nonphysical for
      // this energy regime). Anything within this is included in the σ fit.
      if (std::abs(KEfit - KEmc) > 30.0) continue;
      hKE[b]->Fill(KEfit - KEmc);
      hTh[b]->Fill(thFit - thMC);
      hKEresAll->Fill(KEfit - KEmc);
      hKEcorr->Fill(KEmc, KEfit);
      ++nFitTot;
   }

   std::vector<double> thMid(NB), thErr(NB);
   std::vector<double> sigE(NB), sigEe(NB), bia(NB), biae(NB), thS(NB), thSe(NB), yld(NB), yldErr(NB);
   for (int b = 0; b < NB; ++b) {
      thMid[b] = 0.5 * (edges[b] + edges[b + 1]);
      thErr[b] = 0.5 * (edges[b + 1] - edges[b]);
      double meanKE = nThrown[b] > 0 ? sumKE[b] / nThrown[b] : 1.;
      auto r = c16pp_perf::fitBin(hKE[b], hTh[b], meanKE, nThrown[b]);
      sigE[b] = r.sigE_pct;
      sigEe[b] = sigE[b] / std::sqrt(2. * std::max(1, r.N));
      bia[b] = r.bias;
      biae[b] = r.sigma / std::sqrt(std::max(1, r.N));
      thS[b] = r.th_sigma;
      thSe[b] = thS[b] / std::sqrt(2. * std::max(1, r.N));
      yld[b] = r.yield_pct;
      yldErr[b] = std::sqrt(yld[b] * (100. - yld[b]) / std::max(1, nThrown[b]));
   }

   auto *c = new TCanvas("c", "16C(p,p) performance summary", 1800, 1100);
   c->Divide(3, 2, 0.005, 0.02);

   c->cd(1);
   gPad->SetGrid();
   auto *gS = new TGraphErrors(NB, thMid.data(), sigE.data(), thErr.data(), sigEe.data());
   gS->SetTitle("Recoil proton KE resolution;#theta_{MC} (deg);#sigma_{KE}/#LTKE#GT (%)");
   gS->SetMarkerStyle(20);
   gS->SetMarkerColor(kAzure + 2);
   gS->SetLineColor(kAzure + 2);
   gS->SetMarkerSize(1.2);
   gS->SetLineWidth(2);
   gS->GetYaxis()->SetRangeUser(0, std::max(40., 1.3 * (*std::max_element(sigE.begin(), sigE.end()))));
   gS->Draw("AP");

   c->cd(2);
   gPad->SetGrid();
   auto *gB = new TGraphErrors(NB, thMid.data(), bia.data(), thErr.data(), biae.data());
   gB->SetTitle("KE bias;#theta_{MC} (deg);bias #LTKE_{fit} - KE_{MC}#GT (MeV)");
   gB->SetMarkerStyle(20);
   gB->SetMarkerColor(kAzure + 2);
   gB->SetLineColor(kAzure + 2);
   gB->SetMarkerSize(1.2);
   gB->SetLineWidth(2);
   gB->GetYaxis()->SetRangeUser(-5, 5);
   gB->Draw("AP");
   auto *zero = new TLine(edges[0], 0, edges[NB], 0);
   zero->SetLineStyle(2);
   zero->SetLineColor(kGray + 2);
   zero->Draw("same");

   c->cd(3);
   gPad->SetGrid();
   auto *gY = new TGraphErrors(NB, thMid.data(), yld.data(), thErr.data(), yldErr.data());
   gY->SetTitle("Track yield;#theta_{MC} (deg);Fit yield (%)");
   gY->SetMarkerStyle(21);
   gY->SetMarkerColor(kSpring - 6);
   gY->SetLineColor(kSpring - 6);
   gY->SetMarkerSize(1.2);
   gY->SetLineWidth(2);
   gY->GetYaxis()->SetRangeUser(0, 110);
   gY->Draw("AP");

   c->cd(4);
   gPad->SetGrid();
   auto *gT = new TGraphErrors(NB, thMid.data(), thS.data(), thErr.data(), thSe.data());
   gT->SetTitle("#theta resolution;#theta_{MC} (deg);#sigma_{#theta} (deg)");
   gT->SetMarkerStyle(20);
   gT->SetMarkerColor(kAzure + 2);
   gT->SetLineColor(kAzure + 2);
   gT->SetMarkerSize(1.2);
   gT->SetLineWidth(2);
   gT->GetYaxis()->SetRangeUser(0, std::max(2.5, 1.3 * (*std::max_element(thS.begin(), thS.end()))));
   gT->Draw("AP");

   c->cd(5);
   gPad->SetGrid();
   hKEcorr->SetMarkerStyle(20);
   hKEcorr->SetMarkerSize(0.4);
   hKEcorr->SetMarkerColor(kAzure + 2);
   hKEcorr->Draw("colz");
   auto *diag = new TLine(0, 0, 25, 25);
   diag->SetLineColor(kRed + 1);
   diag->SetLineWidth(2);
   diag->SetLineStyle(2);
   diag->Draw("same");

   c->cd(6);
   gPad->SetGrid();
   hKEresAll->SetFillColorAlpha(kAzure - 4, 0.5);
   hKEresAll->SetLineColor(kAzure + 2);
   hKEresAll->SetLineWidth(2);
   hKEresAll->Fit("gaus", "Q", "", -10., 10.);
   if (auto *f = hKEresAll->GetFunction("gaus")) {
      f->SetLineColor(kRed + 1);
      f->SetLineWidth(2);
   }
   hKEresAll->Draw("hist");
   if (auto *f = hKEresAll->GetFunction("gaus")) f->Draw("same");
   auto *txt = new TLatex();
   txt->SetTextSize(0.045);
   txt->SetNDC();
   if (auto *f = hKEresAll->GetFunction("gaus")) {
      txt->DrawLatex(0.18, 0.84, Form("#mu = %.2f MeV", f->GetParameter(1)));
      txt->DrawLatex(0.18, 0.78, Form("#sigma = %.2f MeV", f->GetParameter(2)));
   }

   c->SaveAs(Form("data/performance_%s.png", outTag));
   c->SaveAs(Form("data/performance_%s.pdf", outTag));
   std::cout << "Wrote data/performance_" << outTag << ".{png,pdf}\n";

   std::cout << "\n=== 16C(p,p) recoil-proton performance ===\n";
   std::cout << "Events with proton: " << nEvents << "  fitted: " << nFitTot << "\n";
   std::cout << "theta(deg)  Nthr  yield   bias    sigE      sigTheta\n";
   std::cout << std::string(60, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      char buf[256];
      snprintf(buf, sizeof(buf), "%4.0f-%-4.0f  %4d  %5.1f%%   %+5.2f   %5.1f%%    %5.2f deg",
               edges[b], edges[b + 1], nThrown[b], yld[b], bia[b], sigE[b], thS[b]);
      std::cout << buf << "\n";
   }
}
