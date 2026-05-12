/// @file make_performance.C
/// @brief One-shot performance summary for the pi-TPC sandbox.
///
/// Reads ./data/attpcsim.root, ./data/output_ukf_only.root (PRA+UKF), and
/// optionally ./data/output_ukf_truthpra.root (intrinsic UKF ceiling), and
/// produces ./data/performance.png/pdf — a 6-panel summary suitable for
/// sharing.
///
/// Run from pi_TPC/: `root -b -q analysis/make_performance.C`

namespace pi_perf {
struct BinResult {
   int N;
   double bias, sigma, sigE_pct;
   double th_sigma, yield_pct;
};

BinResult fitBin(TH1F *hKE, TH1F *hTh, double meanKE, int nThrown)
{
   BinResult r{0, 0, 0, 0, 0, 0};
   if (hKE->GetEntries() < 15) return r;
   double rms = hKE->GetRMS(), mean = hKE->GetMean();
   double a = std::max(-25., mean - 2.5 * rms), b = std::min(25., mean + 2.5 * rms);
   hKE->Fit("gaus", "Q0", "", a, b);
   auto *f = hKE->GetFunction("gaus");
   r.bias = f ? f->GetParameter(1) : mean;
   r.sigma = f ? f->GetParameter(2) : rms;
   r.sigE_pct = (meanKE > 0 && std::isfinite(r.sigma)) ? r.sigma / meanKE * 100. : 0;
   if (hTh->GetEntries() >= 15) {
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
} // namespace pi_perf

void make_performance(const char *ukfFile = "data/output_ukf_only.root",
                      const char *idealFile = "data/output_ukf_truthpra.root",
                      const char *outPrefix = "data/performance",
                      const char *titleTag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim("data/attpcsim.root");
   TFile fPRA(ukfFile);
   TFile fIDL(idealFile, "READ");
   bool haveIdeal = !fIDL.IsZombie() && fIDL.Get("cbmsim") != nullptr;

   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tPRA = (TTree *)fPRA.Get("cbmsim");
   auto *tIDL = haveIdeal ? (TTree *)fIDL.Get("cbmsim") : nullptr;

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teP = new TClonesArray("AtTrackingEvent");
   TClonesArray *teI = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tPRA->SetBranchAddress("AtTrackingEvent", &teP);
   if (tIDL) tIDL->SetBranchAddress("AtTrackingEvent", &teI);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};

   std::vector<TH1F *> hKEp(NB), hThp(NB), hKEi(NB), hThi(NB);
   std::vector<int> nThrown(NB, 0);
   std::vector<double> sumKE(NB, 0.);
   for (int b = 0; b < NB; ++b) {
      hKEp[b] = new TH1F(Form("hKEp_%d", b), "", 80, -25., 25.);
      hThp[b] = new TH1F(Form("hThp_%d", b), "", 80, -10., 10.);
      hKEi[b] = new TH1F(Form("hKEi_%d", b), "", 80, -25., 25.);
      hThi[b] = new TH1F(Form("hThi_%d", b), "", 80, -10., 10.);
   }

   // Global histograms
   TH1F *hKEresAll = new TH1F("hKEresAll", "KE residual (all #theta);KE_{fit} - KE_{MC} (MeV);counts", 80, -25., 25.);
   TH2F *hKEcorr = new TH2F("hKEcorr", "KE_{fit} vs KE_{MC};KE_{MC} (MeV);KE_{fit} (MeV)", 60, 0., 55., 60, 0., 55.);
   TH1F *hYieldThr = new TH1F("hYieldThr", ";#theta_{MC} (deg);counts", NB, edges);
   TH1F *hYieldFit = new TH1F("hYieldFit", ";#theta_{MC} (deg);counts", NB, edges);

   auto fillFromTree = [&](TTree *t, TClonesArray *te, std::vector<TH1F *> &hKE, std::vector<TH1F *> &hTh,
                            bool isPRA) {
      Long64_t n = std::min(tSim->GetEntries(), t->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         tSim->GetEntry(i);
         t->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;
         double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
         double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
         double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
         int b = -1;
         for (int k = 0; k < NB; ++k) if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
         if (b < 0) continue;
         if (isPRA) { ++nThrown[b]; sumKE[b] += KEmc; hYieldThr->Fill(thMC); }

         if (te->GetEntries() == 0) continue;
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
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
         hKE[b]->Fill(KEfit - KEmc);
         hTh[b]->Fill(thFit - thMC);
         if (isPRA) {
            hKEresAll->Fill(KEfit - KEmc);
            hKEcorr->Fill(KEmc, KEfit);
            hYieldFit->Fill(thMC);
         }
      }
   };
   fillFromTree(tPRA, teP, hKEp, hThp, true);
   if (tIDL) fillFromTree(tIDL, teI, hKEi, hThi, false);

   std::vector<double> thMid(NB), thErr(NB);
   std::vector<double> sigP(NB), sigPe(NB), sigI(NB), sigIe(NB);
   std::vector<double> biaP(NB), biaPe(NB), biaI(NB), biaIe(NB);
   std::vector<double> thsP(NB), thsPe(NB), thsI(NB), thsIe(NB);
   std::vector<double> yld(NB), yldErr(NB);
   for (int b = 0; b < NB; ++b) {
      thMid[b] = 0.5 * (edges[b] + edges[b + 1]);
      thErr[b] = 0.5 * (edges[b + 1] - edges[b]);
      double meanKE = nThrown[b] > 0 ? sumKE[b] / nThrown[b] : 1.;
      auto rP = pi_perf::fitBin(hKEp[b], hThp[b], meanKE, nThrown[b]);
      sigP[b] = rP.sigE_pct;
      sigPe[b] = sigP[b] / std::sqrt(2. * std::max(1, rP.N));
      biaP[b] = rP.bias;
      biaPe[b] = rP.sigma / std::sqrt(std::max(1, rP.N));
      thsP[b] = rP.th_sigma;
      thsPe[b] = thsP[b] / std::sqrt(2. * std::max(1, rP.N));
      yld[b] = rP.yield_pct;
      yldErr[b] = std::sqrt(yld[b] * (100. - yld[b]) / std::max(1, nThrown[b]));
      if (tIDL) {
         auto rI = pi_perf::fitBin(hKEi[b], hThi[b], meanKE, nThrown[b]);
         sigI[b] = rI.sigE_pct;
         sigIe[b] = sigI[b] / std::sqrt(2. * std::max(1, rI.N));
         biaI[b] = rI.bias;
         biaIe[b] = rI.sigma / std::sqrt(std::max(1, rI.N));
         thsI[b] = rI.th_sigma;
         thsIe[b] = thsI[b] / std::sqrt(2. * std::max(1, rI.N));
      }
   }

   auto *c = new TCanvas("c", "pi-TPC performance summary", 1800, 1100);
   c->Divide(3, 2, 0.005, 0.02);

   // ---- Panel 1: σ/E vs θ ----
   c->cd(1);
   gPad->SetGrid();
   auto *gSp = new TGraphErrors(NB, thMid.data(), sigP.data(), thErr.data(), sigPe.data());
   gSp->SetTitle("KE resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{KE}/#LTKE#GT (%)");
   gSp->SetMarkerStyle(20);
   gSp->SetMarkerColor(kAzure + 2);
   gSp->SetMarkerSize(1.2);
   gSp->SetLineColor(kAzure + 2);
   gSp->SetLineWidth(2);
   gSp->GetYaxis()->SetRangeUser(0, 35);
   gSp->Draw("AP");
   if (tIDL) {
      auto *gSi = new TGraphErrors(NB, thMid.data(), sigI.data(), thErr.data(), sigIe.data());
      gSi->SetMarkerStyle(24);
      gSi->SetMarkerColor(kRed + 1);
      gSi->SetMarkerSize(1.2);
      gSi->SetLineColor(kRed + 1);
      gSi->SetLineWidth(2);
      gSi->Draw("P same");
      auto *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(gSp, "PRA + UKF", "lp");
      leg->AddEntry(gSi, "Truth seed (intrinsic)", "lp");
      leg->Draw();
   }

   // ---- Panel 2: KE bias vs θ ----
   c->cd(2);
   gPad->SetGrid();
   auto *gBp = new TGraphErrors(NB, thMid.data(), biaP.data(), thErr.data(), biaPe.data());
   gBp->SetTitle("KE bias vs #theta_{MC};#theta_{MC} (deg);bias #LTKE_{fit} - KE_{MC}#GT (MeV)");
   gBp->SetMarkerStyle(20);
   gBp->SetMarkerColor(kAzure + 2);
   gBp->SetMarkerSize(1.2);
   gBp->SetLineColor(kAzure + 2);
   gBp->SetLineWidth(2);
   gBp->GetYaxis()->SetRangeUser(-5, 5);
   gBp->Draw("AP");
   if (tIDL) {
      auto *gBi = new TGraphErrors(NB, thMid.data(), biaI.data(), thErr.data(), biaIe.data());
      gBi->SetMarkerStyle(24);
      gBi->SetMarkerColor(kRed + 1);
      gBi->SetMarkerSize(1.2);
      gBi->SetLineColor(kRed + 1);
      gBi->SetLineWidth(2);
      gBi->Draw("P same");
   }
   TLine *zero = new TLine(edges[0], 0, edges[NB], 0);
   zero->SetLineStyle(2);
   zero->SetLineColor(kGray + 2);
   zero->Draw("same");

   // ---- Panel 3: PRA+UKF yield vs θ ----
   c->cd(3);
   gPad->SetGrid();
   auto *gY = new TGraphErrors(NB, thMid.data(), yld.data(), thErr.data(), yldErr.data());
   gY->SetTitle("Track yield vs #theta_{MC};#theta_{MC} (deg);Fit yield (%)");
   gY->SetMarkerStyle(21);
   gY->SetMarkerColor(kSpring - 6);
   gY->SetMarkerSize(1.2);
   gY->SetLineColor(kSpring - 6);
   gY->SetLineWidth(2);
   gY->GetYaxis()->SetRangeUser(0, 110);
   gY->Draw("AP");

   // ---- Panel 4: θ resolution vs θ ----
   c->cd(4);
   gPad->SetGrid();
   auto *gTp = new TGraphErrors(NB, thMid.data(), thsP.data(), thErr.data(), thsPe.data());
   gTp->SetTitle("#theta resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{#theta} (deg)");
   gTp->SetMarkerStyle(20);
   gTp->SetMarkerColor(kAzure + 2);
   gTp->SetMarkerSize(1.2);
   gTp->SetLineColor(kAzure + 2);
   gTp->SetLineWidth(2);
   gTp->GetYaxis()->SetRangeUser(0, 1.5);
   gTp->Draw("AP");
   if (tIDL) {
      auto *gTi = new TGraphErrors(NB, thMid.data(), thsI.data(), thErr.data(), thsIe.data());
      gTi->SetMarkerStyle(24);
      gTi->SetMarkerColor(kRed + 1);
      gTi->SetMarkerSize(1.2);
      gTi->SetLineColor(kRed + 1);
      gTi->SetLineWidth(2);
      gTi->Draw("P same");
   }

   // ---- Panel 5: KE_fit vs KE_MC scatter ----
   c->cd(5);
   gPad->SetGrid();
   hKEcorr->SetMarkerStyle(20);
   hKEcorr->SetMarkerSize(0.4);
   hKEcorr->SetMarkerColor(kAzure + 2);
   hKEcorr->Draw("colz");
   auto *diag = new TLine(0, 0, 55, 55);
   diag->SetLineColor(kRed + 1);
   diag->SetLineWidth(2);
   diag->SetLineStyle(2);
   diag->Draw("same");

   // ---- Panel 6: global KE residual + Gaussian core fit ----
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

   TString png = TString(outPrefix) + ".png";
   TString pdf = TString(outPrefix) + ".pdf";
   c->SaveAs(png);
   c->SaveAs(pdf);
   if (std::strlen(titleTag) > 0)
      c->SetTitle(titleTag);
   std::cout << "Wrote " << outPrefix << ".{png,pdf}\n";

   // Print summary table
   std::cout << "\n=== pi-TPC performance summary ===\n";
   std::cout << "theta(deg)  Nthr  yield   bias_PRA  sigE_PRA";
   if (tIDL) std::cout << "   sigE_truth";
   std::cout << "   sigTheta_PRA\n";
   std::cout << std::string(80, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      char buf[256];
      snprintf(buf, sizeof(buf), "%4.0f-%-4.0f  %4d  %5.1f%%   %+5.2f      %5.1f%%",
               edges[b], edges[b + 1], nThrown[b], yld[b], biaP[b], sigP[b]);
      std::cout << buf;
      if (tIDL) printf("    %5.1f%%", sigI[b]);
      printf("        %5.2f deg\n", thsP[b]);
   }
}
