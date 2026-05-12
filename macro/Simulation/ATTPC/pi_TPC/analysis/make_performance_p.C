/// @file make_performance_p.C
/// @brief 6-panel performance summary for relativistic pions where σ_p/p
/// (not σ_KE/KE) is the right figure of merit. Mirrors make_performance.C
/// but plots momentum residual instead of KE residual.
///
/// Usage:
///   root -b -q 'analysis/make_performance_p.C("data/output_ukf_only_800MeV.root",
///       "data/output_ukf_truthpra_800MeV.root","data/attpcsim_800MeV.root",
///       "data/results_800MeV_2T_2026-05-13/performance_p","800 MeV/c, B=2 T")'

namespace pi_perf_p {
struct BinResult {
   int N;
   double bias, sigma, sigP_pct;
   double th_sigma, yield_pct;
};

BinResult fitBin(TH1F *hDP, TH1F *hTh, int nThrown)
{
   BinResult r{0, 0, 0, 0, 0, 0};
   if (hDP->GetEntries() < 15) return r;
   double rms = hDP->GetRMS(), mean = hDP->GetMean();
   hDP->Fit("gaus", "Q0", "", mean - 2.0 * rms, mean + 2.0 * rms);
   auto *f = hDP->GetFunction("gaus");
   r.bias = f ? f->GetParameter(1) : mean;
   r.sigma = f ? f->GetParameter(2) : rms;
   r.sigP_pct = r.sigma * 100.;
   if (hTh->GetEntries() >= 15) {
      double rT = hTh->GetRMS(), mT = hTh->GetMean();
      hTh->Fit("gaus", "Q0", "", mT - 2.5 * rT, mT + 2.5 * rT);
      auto *fT = hTh->GetFunction("gaus");
      r.th_sigma = fT ? fT->GetParameter(2) : rT;
   }
   r.N = hDP->GetEntries();
   r.yield_pct = nThrown > 0 ? 100. * r.N / nThrown : 0;
   return r;
}
} // namespace

void make_performance_p(const char *ukfFile = "data/output_ukf_only_800MeV.root",
                        const char *idealFile = "data/output_ukf_truthpra_800MeV.root",
                        const char *simFile = "data/attpcsim_800MeV.root",
                        const char *outPrefix = "data/performance_p",
                        const char *titleTag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim(simFile);
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

   std::vector<TH1F *> hDPp(NB), hThp(NB), hDPi(NB), hThi(NB);
   std::vector<int> nThrown(NB, 0);
   for (int b = 0; b < NB; ++b) {
      hDPp[b] = new TH1F(Form("hDPp_%d", b), "", 80, -0.5, 0.5);
      hThp[b] = new TH1F(Form("hThp_%d", b), "", 80, -3., 3.);
      hDPi[b] = new TH1F(Form("hDPi_%d", b), "", 80, -0.1, 0.1);
      hThi[b] = new TH1F(Form("hThi_%d", b), "", 80, -3., 3.);
   }

   TH1F *hDPresAll = new TH1F("hDPresAll", "p residual (all #theta);(p_{fit}-p_{MC})/p_{MC};counts", 80, -0.6, 0.6);
   TH2F *hPcorr = new TH2F("hPcorr", "p_{fit} vs p_{MC};p_{MC} (MeV/c);p_{fit} (MeV/c)", 60, 500., 900., 60, 500., 900.);
   TH1F *hYieldThr = new TH1F("hYieldThr", ";#theta_{MC} (deg);counts", NB, edges);
   TH1F *hYieldFit = new TH1F("hYieldFit", ";#theta_{MC} (deg);counts", NB, edges);

   const double mass_pi = 139.57039;

   auto fillFromTree = [&](TTree *t, TClonesArray *te, std::vector<TH1F *> &hDP, std::vector<TH1F *> &hTh, bool isPRA) {
      Long64_t n = std::min(tSim->GetEntries(), t->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         tSim->GetEntry(i);
         t->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;
         double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz()) * 1000.;
         double thMC = std::acos(mc->GetPz() / (pmc / 1000.)) * 180. / M_PI;
         int b = -1;
         for (int k = 0; k < NB; ++k)
            if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
         if (b < 0) continue;
         if (isPRA) { ++nThrown[b]; hYieldThr->Fill(thMC); }
         if (te->GetEntries() == 0) continue;
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &tr : fitted) {
            if (!tr->GetTrackMetadata()) continue;
            double ndf = tr->GetTrackMetadata()->GetNdf();
            double c = ndf > 0 ? tr->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (c < bestChi) { bestChi = c; best = tr.get(); }
         }
         if (!best) continue;
         double KEfit = best->GetKinematics().kineticEnergy;
         double Efit = KEfit + mass_pi;
         double pfit = std::sqrt(Efit * Efit - mass_pi * mass_pi);
         double thFit = best->GetKinematics().theta * 180. / M_PI;
         hDP[b]->Fill((pfit - pmc) / pmc);
         hTh[b]->Fill(thFit - thMC);
         if (isPRA) {
            hDPresAll->Fill((pfit - pmc) / pmc);
            hPcorr->Fill(pmc, pfit);
            hYieldFit->Fill(thMC);
         }
      }
   };
   fillFromTree(tPRA, teP, hDPp, hThp, true);
   if (tIDL) fillFromTree(tIDL, teI, hDPi, hThi, false);

   std::vector<double> thMid(NB), thErr(NB);
   std::vector<double> sigP(NB), sigPe(NB), sigI(NB), sigIe(NB);
   std::vector<double> biaP(NB), biaPe(NB), biaI(NB), biaIe(NB);
   std::vector<double> thsP(NB), thsPe(NB), thsI(NB), thsIe(NB);
   std::vector<double> yld(NB), yldErr(NB);
   for (int b = 0; b < NB; ++b) {
      thMid[b] = 0.5 * (edges[b] + edges[b + 1]);
      thErr[b] = 0.5 * (edges[b + 1] - edges[b]);
      auto rP = pi_perf_p::fitBin(hDPp[b], hThp[b], nThrown[b]);
      sigP[b] = rP.sigP_pct;
      sigPe[b] = sigP[b] / std::sqrt(2. * std::max(1, rP.N));
      biaP[b] = rP.bias * 100;
      biaPe[b] = rP.sigma * 100 / std::sqrt(std::max(1, rP.N));
      thsP[b] = rP.th_sigma;
      thsPe[b] = thsP[b] / std::sqrt(2. * std::max(1, rP.N));
      yld[b] = rP.yield_pct;
      yldErr[b] = std::sqrt(yld[b] * (100. - yld[b]) / std::max(1, nThrown[b]));
      if (tIDL) {
         auto rI = pi_perf_p::fitBin(hDPi[b], hThi[b], nThrown[b]);
         sigI[b] = rI.sigP_pct;
         sigIe[b] = sigI[b] / std::sqrt(2. * std::max(1, rI.N));
         biaI[b] = rI.bias * 100;
         biaIe[b] = rI.sigma * 100 / std::sqrt(std::max(1, rI.N));
         thsI[b] = rI.th_sigma;
         thsIe[b] = thsI[b] / std::sqrt(2. * std::max(1, rI.N));
      }
   }

   auto *c = new TCanvas("c", titleTag, 1800, 1100);
   c->Divide(3, 2, 0.005, 0.02);

   c->cd(1); gPad->SetGrid();
   auto *gSp = new TGraphErrors(NB, thMid.data(), sigP.data(), thErr.data(), sigPe.data());
   gSp->SetTitle("p resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{p}/p (%)");
   gSp->SetMarkerStyle(20); gSp->SetMarkerColor(kAzure + 2); gSp->SetLineColor(kAzure + 2);
   gSp->SetMarkerSize(1.2); gSp->SetLineWidth(2);
   gSp->GetYaxis()->SetRangeUser(0, 20);
   gSp->Draw("AP");
   if (tIDL) {
      auto *gSi = new TGraphErrors(NB, thMid.data(), sigI.data(), thErr.data(), sigIe.data());
      gSi->SetMarkerStyle(24); gSi->SetMarkerColor(kRed + 1); gSi->SetLineColor(kRed + 1);
      gSi->SetMarkerSize(1.2); gSi->SetLineWidth(2);
      gSi->Draw("P same");
      auto *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry(gSp, "PRA + UKF", "lp");
      leg->AddEntry(gSi, "Truth seed (intrinsic)", "lp");
      leg->Draw();
   }

   c->cd(2); gPad->SetGrid();
   auto *gBp = new TGraphErrors(NB, thMid.data(), biaP.data(), thErr.data(), biaPe.data());
   gBp->SetTitle("p bias vs #theta_{MC};#theta_{MC} (deg);bias #LT(p_{fit}-p_{MC})/p#GT (%)");
   gBp->SetMarkerStyle(20); gBp->SetMarkerColor(kAzure + 2); gBp->SetLineColor(kAzure + 2);
   gBp->SetMarkerSize(1.2); gBp->SetLineWidth(2);
   gBp->GetYaxis()->SetRangeUser(-5, 5);
   gBp->Draw("AP");
   if (tIDL) {
      auto *gBi = new TGraphErrors(NB, thMid.data(), biaI.data(), thErr.data(), biaIe.data());
      gBi->SetMarkerStyle(24); gBi->SetMarkerColor(kRed + 1); gBi->SetLineColor(kRed + 1);
      gBi->SetMarkerSize(1.2); gBi->SetLineWidth(2);
      gBi->Draw("P same");
   }
   TLine *zero = new TLine(edges[0], 0, edges[NB], 0);
   zero->SetLineStyle(2); zero->SetLineColor(kGray + 2); zero->Draw("same");

   c->cd(3); gPad->SetGrid();
   auto *gY = new TGraphErrors(NB, thMid.data(), yld.data(), thErr.data(), yldErr.data());
   gY->SetTitle("Track yield vs #theta_{MC};#theta_{MC} (deg);Fit yield (%)");
   gY->SetMarkerStyle(21); gY->SetMarkerColor(kSpring - 6); gY->SetLineColor(kSpring - 6);
   gY->SetMarkerSize(1.2); gY->SetLineWidth(2);
   gY->GetYaxis()->SetRangeUser(0, 110);
   gY->Draw("AP");

   c->cd(4); gPad->SetGrid();
   auto *gTp = new TGraphErrors(NB, thMid.data(), thsP.data(), thErr.data(), thsPe.data());
   gTp->SetTitle("#theta resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{#theta} (deg)");
   gTp->SetMarkerStyle(20); gTp->SetMarkerColor(kAzure + 2); gTp->SetLineColor(kAzure + 2);
   gTp->SetMarkerSize(1.2); gTp->SetLineWidth(2);
   gTp->GetYaxis()->SetRangeUser(0, 1.5);
   gTp->Draw("AP");
   if (tIDL) {
      auto *gTi = new TGraphErrors(NB, thMid.data(), thsI.data(), thErr.data(), thsIe.data());
      gTi->SetMarkerStyle(24); gTi->SetMarkerColor(kRed + 1); gTi->SetLineColor(kRed + 1);
      gTi->SetMarkerSize(1.2); gTi->SetLineWidth(2);
      gTi->Draw("P same");
   }

   c->cd(5); gPad->SetGrid();
   hPcorr->SetMarkerStyle(20); hPcorr->SetMarkerSize(0.4); hPcorr->SetMarkerColor(kAzure + 2);
   hPcorr->Draw("colz");
   auto *diag = new TLine(500, 500, 900, 900);
   diag->SetLineColor(kRed + 1); diag->SetLineWidth(2); diag->SetLineStyle(2);
   diag->Draw("same");

   c->cd(6); gPad->SetGrid();
   hDPresAll->SetFillColorAlpha(kAzure - 4, 0.5);
   hDPresAll->SetLineColor(kAzure + 2);
   hDPresAll->SetLineWidth(2);
   hDPresAll->Fit("gaus", "Q", "", -0.2, 0.2);
   if (auto *f = hDPresAll->GetFunction("gaus")) {
      f->SetLineColor(kRed + 1); f->SetLineWidth(2);
   }
   hDPresAll->Draw("hist");
   if (auto *f = hDPresAll->GetFunction("gaus")) f->Draw("same");
   auto *txt = new TLatex();
   txt->SetTextSize(0.045); txt->SetNDC();
   if (auto *f = hDPresAll->GetFunction("gaus")) {
      txt->DrawLatex(0.18, 0.84, Form("#mu = %.2f%%", f->GetParameter(1) * 100));
      txt->DrawLatex(0.18, 0.78, Form("#sigma = %.2f%%", f->GetParameter(2) * 100));
   }

   TString png = TString(outPrefix) + ".png";
   TString pdf = TString(outPrefix) + ".pdf";
   c->SaveAs(png);
   c->SaveAs(pdf);
   std::cout << "Wrote " << outPrefix << ".{png,pdf}\n";

   std::cout << "\n=== σ_p/p summary [" << titleTag << "] ===\n";
   std::cout << "theta(deg)  Nthr  yield   bias_PRA  σ_p/p_PRA";
   if (tIDL) std::cout << "   σ_p/p_truth";
   std::cout << "   σ_θ_PRA\n";
   std::cout << std::string(80, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      printf("%4.0f-%-4.0f  %4d  %5.1f%%   %+5.2f%%    %5.2f%%", edges[b], edges[b + 1], nThrown[b], yld[b], biaP[b], sigP[b]);
      if (tIDL) printf("    %5.2f%%", sigI[b]);
      printf("        %5.2f deg\n", thsP[b]);
   }
}
