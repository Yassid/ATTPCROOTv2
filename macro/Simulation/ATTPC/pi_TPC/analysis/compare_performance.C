/// @file compare_performance.C
/// @brief Overlay σ/E and bias vs θ for two UKF outputs (e.g. baseline vs
/// arc-aligned pre-cluster). Reads MC truth from attpcsim.root, applies the
/// same fitBin logic as analysis/make_performance.C.
///
/// Usage:
///   root -b -q 'analysis/compare_performance.C("data/output_ukf_only.root",
///       "data/output_ukf_only_arc6.root","baseline","arc-bin 6 mm",
///       "data/results_arc6_2026-05-13/compare")'

namespace pi_cmp {
struct BinResult { int N; double bias, sigma, sigE_pct, th_sigma; };
BinResult fitBin(TH1F *hKE, TH1F *hTh, double meanKE) {
   BinResult r{0, 0, 0, 0, 0};
   if (hKE->GetEntries() < 15) return r;
   double rms = hKE->GetRMS(), mean = hKE->GetMean();
   hKE->Fit("gaus", "Q0", "", mean - 2.5 * rms, mean + 2.5 * rms);
   auto *f = hKE->GetFunction("gaus");
   r.bias = f ? f->GetParameter(1) : mean;
   r.sigma = f ? f->GetParameter(2) : rms;
   r.sigE_pct = (meanKE > 0 && std::isfinite(r.sigma)) ? r.sigma / meanKE * 100. : 0;
   if (hTh->GetEntries() >= 15) {
      double rT = hTh->GetRMS(), mT = hTh->GetMean();
      hTh->Fit("gaus", "Q0", "", mT - 2.5 * rT, mT + 2.5 * rT);
      auto *fT = hTh->GetFunction("gaus");
      r.th_sigma = fT ? fT->GetParameter(2) : rT;
   }
   r.N = hKE->GetEntries();
   return r;
}

void fill(TTree *tSim, TTree *tUKF, TClonesArray *trks, TClonesArray *te,
          const std::vector<TH1F *> &hKE, const std::vector<TH1F *> &hTh,
          std::vector<double> &sumKE, std::vector<int> &nThr,
          const double *edges, int NB) {
   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
      double pmc = std::sqrt(mc->GetPx()*mc->GetPx()+mc->GetPy()*mc->GetPy()+mc->GetPz()*mc->GetPz());
      double th = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      int b = -1;
      for (int k = 0; k < NB; ++k) if (th >= edges[k] && th < edges[k+1]) { b = k; break; }
      if (b < 0) continue;
      ++nThr[b]; sumKE[b] += KEmc;
      if (te->GetEntries() == 0) continue;
      auto *trkEvt = (AtTrackingEvent *)te->At(0);
      auto &fitted = trkEvt->GetFittedTracks();
      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double c = ndf > 0 ? t->GetTrackMetadata()->GetChi2()/ndf : 1e30;
         if (c < bestChi) { bestChi = c; best = t.get(); }
      }
      if (!best) continue;
      hKE[b]->Fill(best->GetKinematics().kineticEnergy - KEmc);
      hTh[b]->Fill(best->GetKinematics().theta * 180./M_PI - th);
   }
}
} // namespace

void compare_performance(const char *fileA = "data/output_ukf_only.root",
                         const char *fileB = "data/output_ukf_only_arc6.root",
                         const char *labelA = "baseline",
                         const char *labelB = "arc-bin 6 mm",
                         const char *outPrefix = "data/compare")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim("data/attpcsim.root");
   TFile fA(fileA);
   TFile fB(fileB);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tA = (TTree *)fA.Get("cbmsim");
   auto *tB = (TTree *)fB.Get("cbmsim");

   auto *trks = new TClonesArray("AtMCTrack");
   auto *teA = new TClonesArray("AtTrackingEvent");
   auto *teB = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tA->SetBranchAddress("AtTrackingEvent", &teA);
   tB->SetBranchAddress("AtTrackingEvent", &teB);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<TH1F *> hKEa(NB), hTha(NB), hKEb(NB), hThb(NB);
   std::vector<double> sumKEa(NB, 0), sumKEb(NB, 0);
   std::vector<int> nThra(NB, 0), nThrb(NB, 0);
   for (int b = 0; b < NB; ++b) {
      hKEa[b] = new TH1F(Form("hKEa_%d", b), "", 80, -25., 25.);
      hTha[b] = new TH1F(Form("hTha_%d", b), "", 80, -10., 10.);
      hKEb[b] = new TH1F(Form("hKEb_%d", b), "", 80, -25., 25.);
      hThb[b] = new TH1F(Form("hThb_%d", b), "", 80, -10., 10.);
   }
   pi_cmp::fill(tSim, tA, trks, teA, hKEa, hTha, sumKEa, nThra, edges, NB);
   tSim->ResetBranchAddresses();
   tSim->SetBranchAddress("MCTrack", &trks);
   pi_cmp::fill(tSim, tB, trks, teB, hKEb, hThb, sumKEb, nThrb, edges, NB);

   std::vector<double> thMid(NB), thErr(NB);
   std::vector<double> sigA(NB), sigAe(NB), biaA(NB), biaAe(NB), thsA(NB), thsAe(NB);
   std::vector<double> sigB(NB), sigBe(NB), biaB(NB), biaBe(NB), thsB(NB), thsBe(NB);
   for (int b = 0; b < NB; ++b) {
      thMid[b] = 0.5 * (edges[b] + edges[b+1]);
      thErr[b] = 0.5 * (edges[b+1] - edges[b]);
      double meanA = nThra[b] > 0 ? sumKEa[b] / nThra[b] : 1.;
      double meanB = nThrb[b] > 0 ? sumKEb[b] / nThrb[b] : 1.;
      auto rA = pi_cmp::fitBin(hKEa[b], hTha[b], meanA);
      auto rB = pi_cmp::fitBin(hKEb[b], hThb[b], meanB);
      sigA[b]=rA.sigE_pct; sigAe[b]=sigA[b]/std::sqrt(2.*std::max(1,rA.N));
      biaA[b]=rA.bias;     biaAe[b]=rA.sigma/std::sqrt(std::max(1,rA.N));
      thsA[b]=rA.th_sigma; thsAe[b]=thsA[b]/std::sqrt(2.*std::max(1,rA.N));
      sigB[b]=rB.sigE_pct; sigBe[b]=sigB[b]/std::sqrt(2.*std::max(1,rB.N));
      biaB[b]=rB.bias;     biaBe[b]=rB.sigma/std::sqrt(std::max(1,rB.N));
      thsB[b]=rB.th_sigma; thsBe[b]=thsB[b]/std::sqrt(2.*std::max(1,rB.N));
   }

   auto *c = new TCanvas("c", Form("%s vs %s", labelA, labelB), 1500, 500);
   c->Divide(3, 1, 0.005, 0.02);

   // --- Panel 1: sigma/E ---
   c->cd(1);
   gPad->SetGrid();
   auto *gA = new TGraphErrors(NB, thMid.data(), sigA.data(), thErr.data(), sigAe.data());
   auto *gB = new TGraphErrors(NB, thMid.data(), sigB.data(), thErr.data(), sigBe.data());
   gA->SetTitle("KE resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{KE}/#LTKE#GT (%)");
   gA->SetMarkerStyle(20); gA->SetMarkerColor(kAzure+2); gA->SetLineColor(kAzure+2);
   gA->SetMarkerSize(1.2); gA->SetLineWidth(2);
   gA->GetYaxis()->SetRangeUser(0, 35);
   gB->SetMarkerStyle(21); gB->SetMarkerColor(kRed+1); gB->SetLineColor(kRed+1);
   gB->SetMarkerSize(1.2); gB->SetLineWidth(2);
   gA->Draw("AP");
   gB->Draw("P same");
   auto *leg1 = new TLegend(0.40, 0.72, 0.88, 0.88);
   leg1->SetBorderSize(0); leg1->SetFillStyle(0);
   leg1->AddEntry(gA, labelA, "lp"); leg1->AddEntry(gB, labelB, "lp");
   leg1->Draw();

   // --- Panel 2: bias ---
   c->cd(2);
   gPad->SetGrid();
   auto *gBA = new TGraphErrors(NB, thMid.data(), biaA.data(), thErr.data(), biaAe.data());
   auto *gBB = new TGraphErrors(NB, thMid.data(), biaB.data(), thErr.data(), biaBe.data());
   gBA->SetTitle("KE bias vs #theta_{MC};#theta_{MC} (deg);bias #LTKE_{fit}-KE_{MC}#GT (MeV)");
   gBA->SetMarkerStyle(20); gBA->SetMarkerColor(kAzure+2); gBA->SetLineColor(kAzure+2);
   gBA->SetMarkerSize(1.2); gBA->SetLineWidth(2);
   gBA->GetYaxis()->SetRangeUser(-5, 5);
   gBB->SetMarkerStyle(21); gBB->SetMarkerColor(kRed+1); gBB->SetLineColor(kRed+1);
   gBB->SetMarkerSize(1.2); gBB->SetLineWidth(2);
   gBA->Draw("AP");
   gBB->Draw("P same");
   TLine *zero = new TLine(edges[0], 0, edges[NB], 0);
   zero->SetLineStyle(2); zero->SetLineColor(kGray+2); zero->Draw("same");

   // --- Panel 3: theta resolution ---
   c->cd(3);
   gPad->SetGrid();
   auto *gTA = new TGraphErrors(NB, thMid.data(), thsA.data(), thErr.data(), thsAe.data());
   auto *gTB = new TGraphErrors(NB, thMid.data(), thsB.data(), thErr.data(), thsBe.data());
   gTA->SetTitle("#theta resolution vs #theta_{MC};#theta_{MC} (deg);#sigma_{#theta} (deg)");
   gTA->SetMarkerStyle(20); gTA->SetMarkerColor(kAzure+2); gTA->SetLineColor(kAzure+2);
   gTA->SetMarkerSize(1.2); gTA->SetLineWidth(2);
   gTA->GetYaxis()->SetRangeUser(0, 1.5);
   gTB->SetMarkerStyle(21); gTB->SetMarkerColor(kRed+1); gTB->SetLineColor(kRed+1);
   gTB->SetMarkerSize(1.2); gTB->SetLineWidth(2);
   gTA->Draw("AP");
   gTB->Draw("P same");

   TString png = TString(outPrefix) + ".png";
   TString pdf = TString(outPrefix) + ".pdf";
   c->SaveAs(png);
   c->SaveAs(pdf);
   std::cout << "Wrote " << outPrefix << ".{png,pdf}\n";

   std::cout << "\n=== compare: " << labelA << " vs " << labelB << " ===\n";
   std::cout << "theta(deg)   sigE_A   sigE_B    Δ      bias_A    bias_B\n";
   std::cout << std::string(60, '-') << "\n";
   for (int b = 0; b < NB; ++b)
      printf("%4.0f-%-4.0f   %5.1f%%   %5.1f%%   %+5.1f   %+5.2f    %+5.2f\n",
             edges[b], edges[b+1], sigA[b], sigB[b], sigB[b]-sigA[b], biaA[b], biaB[b]);
}
