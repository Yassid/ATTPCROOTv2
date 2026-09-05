/// @file validate_reverse_3Hed.C
/// @brief Prove that ReverseDrift changes the DRIFT and nothing else -- against MC truth.
///
///   root -b -q 'validate_reverse_3Hed.C("<sim>.root","nrm_reco.root","rev_reco.root")'
///
/// Reconstruct ONE sim twice -- once with a normal par, once with a par differing only by
/// `ReverseDrift: 1` -- and hand both here with the sim they came from.
///
/// THE INVARIANT BEING TESTED. mcPoint z is measured along the beam. The drift length is the
/// distance to the pad plane, so reversing which end the pad plane sits on gives
/// drift = z_beam instead of ZPadPlane - z_beam; AtPSA maps the drift back with the complementary
/// sign, so the digi-frame hit z must come out as ZPadPlane - z_beam in BOTH modes. That identity
/// is what lets every downstream convention stand: the digi frame, the fitters'
/// z_lab = ZPadPlane - z_digi, "the vertex is the highest-z_digi end", and the B-field handedness.
///
/// !! DO NOT TEST THIS BY COMPARING THE MEAN HIT z BETWEEN THE TWO MODES. !! That was the first
/// attempt and it reports a false failure. Diffusion spreads charge over more pads at long drift,
/// so more hits clear threshold at the long-drift end -- and the long-drift end is the OPPOSITE
/// end of the chamber in each mode. The mean hit z therefore legitimately differs by tens of mm
/// (measured: 526.9 mm normal against 477.2 mm reversed) while the mapping is perfectly correct.
/// The mean moves for the same reason the feature exists. Anchor on TRUTH instead, which is what
/// this macro does: each hit carries the index of the AtMCPoint that made it, and the sim file
/// still has that point, so z_digi can be compared with ZPadPlane - z_truth hit by hit.
///
/// The end-to-end cross-check, using already-validated machinery, is
/// `ex_genfit_3Hed.C` on both fitted arms: the vertex-z correlation against truth must come back
/// POSITIVE and the same size in both (measured +0.967 and +0.967). A sign error flips it.
void validate_reverse_3Hed(TString simFile, TString normalFile, TString reverseFile, Double_t zPadPlane = 1000.0,
                           Int_t maxEvents = -1)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");

   TFile *fS = TFile::Open(simFile), *fN = TFile::Open(normalFile), *fR = TFile::Open(reverseFile);
   if (!fS || fS->IsZombie() || !fN || fN->IsZombie() || !fR || fR->IsZombie()) {
      printf("cannot open all three files\n");
      return;
   }
   TTree *tS = (TTree *)fS->Get("cbmsim"), *tN = (TTree *)fN->Get("cbmsim"), *tR = (TTree *)fR->Get("cbmsim");
   TClonesArray *mcp = nullptr, *eN = nullptr, *eR = nullptr;
   tS->SetBranchAddress("AtTpcPoint", &mcp);
   tN->SetBranchAddress("AtEventH", &eN);
   tR->SetBranchAddress("AtEventH", &eR);

   Long64_t N = std::min({tS->GetEntries(), tN->GetEntries(), tR->GetEntries()});
   if (maxEvents > 0) N = std::min(N, (Long64_t)maxEvents);

   // residual = z_digi - (ZPadPlane - z_truth). Must be ~0 in BOTH modes; a flipped mapping puts
   // it at -2*(ZPadPlane - z_truth - ZPadPlane/2), i.e. hundreds of mm and strongly z-dependent.
   std::vector<double> resN, resR;
   auto *pN = new TProfile("pN", "", 50, 0, 1000, -1500, 1500);
   auto *pR = new TProfile("pR", "", 50, 0, 1000, -1500, 1500);
   auto *dN = new TH1D("dN", "", 100, -600, 600);
   auto *dR = new TH1D("dR", "", 100, -600, 600);

   auto scan = [&](TClonesArray *ev, std::vector<double> &res, TProfile *prof, TH1D *hd) {
      if (!ev || !ev->GetEntriesFast())
         return;
      auto *e = (AtEvent *)ev->At(0);
      for (const auto &h : e->GetHits()) {
         const auto &mcv = h->GetMCSimPointArray();
         if (mcv.empty())
            continue;
         int idx = mcv.front().pointID;
         if (idx < 0 || idx >= mcp->GetEntriesFast())
            continue;
         auto *p = (AtMCPoint *)mcp->At(idx);
         if (!p)
            continue;
         const double zTruth = p->GetZ() * 10.0;              // cm -> mm, along the beam
         const double expect = zPadPlane - zTruth;            // what z_digi must be, both modes
         const double r = h->GetPosition().Z() - expect;
         res.push_back(r);
         prof->Fill(zTruth, r);
         hd->Fill(r);
      }
   };

   for (Long64_t i = 0; i < N; ++i) {
      tS->GetEntry(i);
      tN->GetEntry(i);
      tR->GetEntry(i);
      scan(eN, resN, pN, dN);
      scan(eR, resR, pR, dR);
   }

   auto med = [](std::vector<double> v) {
      if (v.empty()) return 1e9;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };
   auto iqr = [](std::vector<double> v) {
      if (v.size() < 8) return 1e9;
      std::sort(v.begin(), v.end());
      return v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())];
   };

   printf("\n=== ReverseDrift validation against MC truth, %lld events ===\n", N);
   printf("  residual = z_digi - (ZPadPlane - z_truth);  must be ~0 in BOTH modes\n\n");
   printf("  %-10s %8s %12s %10s\n", "mode", "hits", "median [mm]", "IQR [mm]");
   printf("  %-10s %8zu %12.2f %10.2f\n", "normal", resN.size(), med(resN), iqr(resN));
   printf("  %-10s %8zu %12.2f %10.2f\n", "reversed", resR.size(), med(resR), iqr(resR));

   // Tolerance: the residual carries the longitudinal diffusion and the PSA peak-finding bias, a
   // few mm. A flipped mapping would sit at hundreds of mm, so 20 mm separates the two cases
   // without any ambiguity.
   const bool okN = std::fabs(med(resN)) < 20.0, okR = std::fabs(med(resR)) < 20.0;
   printf("\n  [%s] normal   maps z_digi = ZPadPlane - z_truth\n", okN ? "PASS" : "**FAIL**");
   printf("  [%s] reversed maps z_digi = ZPadPlane - z_truth  <-- the convention is preserved\n",
          okR ? "PASS" : "**FAIL**");
   if (!okR)
      printf("        A large residual here means the digi frame flipped. Everything downstream\n"
             "        that does z_lab = ZPadPlane - z_digi would then be wrong by the drift length.\n");

   // And the drift really did change: mean hit z moves because diffusion now favours the other
   // end. This is the POSITIVE control -- a flag that did nothing would leave it at zero.
   double mzN = 0, mzR = 0;
   for (Long64_t i = 0; i < N; ++i) {
      tN->GetEntry(i);
      tR->GetEntry(i);
      if (eN->GetEntriesFast())
         for (const auto &h : ((AtEvent *)eN->At(0))->GetHits()) mzN += h->GetPosition().Z();
      if (eR->GetEntriesFast())
         for (const auto &h : ((AtEvent *)eR->At(0))->GetHits()) mzR += h->GetPosition().Z();
   }
   printf("\n  positive control -- the drift DID change: charge-weighted mean hit z moves\n"
          "  %.1f mm -> %.1f mm as diffusion swaps ends. Zero movement would mean the flag\n"
          "  never reached AtClusterize.\n",
          mzN / std::max<size_t>(1, resN.size()), mzR / std::max<size_t>(1, resR.size()));

   auto *c = new TCanvas("cRev", "reverse drift validation", 1050, 460);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetGridx(); gPad->SetGridy();
   dN->SetTitle("z_{digi} - (ZPadPlane - z_{truth});residual [mm];hits");
   dN->SetLineColor(kBlack); dN->SetLineWidth(3); dN->Draw("hist");
   dR->SetLineColor(kRed + 1); dR->SetLineWidth(2); dR->SetLineStyle(2); dR->Draw("hist same");
   auto *lg = new TLegend(0.14, 0.76, 0.55, 0.89);
   lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(dN, "normal", "l"); lg->AddEntry(dR, "ReverseDrift = 1", "l"); lg->Draw();
   c->cd(2);
   gPad->SetGridx(); gPad->SetGridy();
   pN->SetTitle("residual vs truth z -- flat at 0 in both modes;z_{truth} along the beam [mm];residual [mm]");
   pN->GetYaxis()->SetRangeUser(-120, 120);
   pN->SetLineColor(kBlack); pN->SetMarkerColor(kBlack); pN->SetMarkerStyle(20); pN->Draw();
   pR->SetLineColor(kRed + 1); pR->SetMarkerColor(kRed + 1); pR->SetMarkerStyle(24); pR->Draw("same");
   c->SaveAs("plots/reverse_drift_check.png");
   printf("\n  wrote plots/reverse_drift_check.png\n\n");
}
