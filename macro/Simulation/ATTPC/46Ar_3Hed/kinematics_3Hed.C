/// @file kinematics_3Hed.C
/// @brief Generated vs reconstructed kinematics for 46Ar(3He,d)47K -- BEFORE ANY FIT.
///
/// Nothing here is fitted. The reconstructed side is AtSpyralPID's circle through the hit cloud,
/// which is all that exists at this stage of the chain:
///     p [MeV/c]  = 299.792458 * Z * Brho [T.m]      (Z = 1 for the deuteron)
///     theta_lab  = 180 - polar                      (the convention AtSpyralPID reports in)
/// Truth comes from MCTrack: the primary deuteron's momentum and start vertex. The same deuteron
/// mass is used on both sides so the comparison carries no offset of its own.
///
/// RUN AT fMinPoints = 15, not the class default 30. At 30 the deuterons between theta_lab 77 and
/// 104 deg are all rejected, which is the middle of the proposal's angular range -- a kinematics
/// plot made at 30 would show two disconnected arcs and hide exactly the region being asked about.
///
/// THE Ex PANEL IS THE POINT. Q value resolution is what the proposal lives or dies by (it asks
/// for ~350 keV FWHM), and Ex can be formed without a fit: two-body inversion of the measured
/// (theta_lab, T_d) against the beam energy AT THE VERTEX,
///     E_R = E_beam(z_vertex) + M(3He) - E_d,   p_R = p_beam - p_d,   Ex = sqrt(E_R^2 - p_R^2) - M(47K)
/// The beam energy at the vertex is not optional: the beam loses 0.957 MeV/cm in this gas
/// (measured, see Ar46_3Hed_sim.C), so 96 MeV across the chamber -- ignoring it would smear Ex by
/// far more than the resolution being measured. What comes out is the PRE-FIT resolution, i.e. the
/// floor a fitter has to improve on, not the proposal's answer.
///
/// HANDEDNESS IS CHECKED, NOT ASSUMED, AND IT IS REVERSED. The simulation mirrors the drift-z
/// sense in digitisation, so the reconstructed vertex comes back at the far end of the chamber:
/// measured correlation of reconstructed against true vertex z is r = -1.000, i.e. an exact
/// mirror, with z_true + z_reco = 100 cm to within the print-out below. The beam energy at the
/// vertex is therefore taken at zUse = DriftLength - z_reco.
///
/// The mirror constant used is the PHYSICAL drift length (100 cm), not a number fitted to the
/// truth -- otherwise the Ex scale would be calibrated on the answer. The printed mean of
/// (z_true + z_reco) is the check on that choice: it has to come out at the drift length, and if
/// it ever does not, this correction is wrong and the Ex offsets below are meaningless.
///
///   root -b -q 'kinematics_3Hed.C()'


/// Analytic two-body locus for the DEUTERON: (theta_lab, T_d) as theta_cm runs 0-180, for a given
/// 47K excitation and a given beam energy. Same relativistic algebra as the generator, so a
/// disagreement between these curves and the generated cloud would mean one of the two is wrong.
///
/// EACH STATE IS A BAND, NOT A LINE. The beam enters at 598 MeV and leaves at ~502 (0.957 MeV/cm
/// over 100 cm), and the reaction happens anywhere along that path, so a state occupies the strip
/// between its entrance and exit curves. Drawing only one line would make the data look worse
/// resolved than it is.
TGraph *ar46_locus(double Ex, double Tb, int color, int style)
{
   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   auto *g = new TGraph();
   double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
   double s = M_b * M_b + M_t * M_t + 2 * Eb * M_t, W = std::sqrt(s);
   double MR = M_R + Ex;
   if (W < MR + M_e) return g;
   double pf = std::sqrt((s - (MR + M_e) * (MR + M_e)) * (s - (MR - M_e) * (MR - M_e))) / (2 * W);
   double Ecm = std::sqrt(pf * pf + M_e * M_e);
   double beta = pb / (Eb + M_t), gam = 1.0 / std::sqrt(1 - beta * beta);
   int n = 0;
   for (double th = 0; th <= 180.0; th += 0.25) {
      double t = th * TMath::DegToRad();
      double pz = gam * (pf * std::cos(t) + beta * Ecm), pt = pf * std::sin(t);
      double E = gam * (Ecm + beta * pf * std::cos(t));
      double thl = std::atan2(pt, pz) * TMath::RadToDeg();
      if (thl < 45 || thl > 145) continue;
      g->SetPoint(n++, thl, E - M_e);
   }
   g->SetLineColor(color);
   g->SetLineStyle(style);
   g->SetLineWidth(2);
   return g;
}

void kinematics_3Hed(TString dir = "/mnt/f/ar46_3hed",
                     TString tags = "gs_s3001,gs_s3002,360_s3011,360_s3012,2020_s3021,2020_s3022",
                     TString png = "plots/kinematics_3Hed.png", Int_t minPoints = 15,
                     Double_t bField = 2.85, Double_t dThetaMax = 10.0, Double_t driftLength = 100.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   // nuclear masses (MeV) -- atomic masses minus Z electrons, as in the generator
   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0;    // MeV at the chamber entrance, 13 MeV/u
   const double dEdz = 0.957;   // MeV/cm, measured beam deposit in this gas
   const double c = 299.792458;

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0)
      spy.SetMinPoints(minPoints);

   auto *hT = new TH2D("hT", "GENERATED (MC truth);#theta_{lab} [deg];T_{d} [MeV]", 100, 50, 140, 100, 0, 70);
   auto *hR = new TH2D("hR", "RECONSTRUCTED, no fit (Spyral B#rho);#theta_{lab} [deg];T_{d} [MeV]",
                       100, 50, 140, 100, 0, 70);
   auto *hV = new TH2D("hV", "vertex z, reco vs truth;truth z [cm];reconstructed z [cm]", 60, 0, 110, 60, -10, 110);
   auto *hEx = new TH1D("hEx", "pre-fit E_{x}(^{47}K);E_{x} [MeV];tracks", 200, -3, 6);
   auto *hEx0 = new TH1D("hEx0", "gs", 200, -3, 6);
   auto *hEx1 = new TH1D("hEx1", "0.36", 200, -3, 6);
   auto *hEx2 = new TH1D("hEx2", "2.02", 200, -3, 6);
   hEx0->SetLineColor(kBlack); hEx1->SetLineColor(kAzure + 2); hEx2->SetLineColor(kRed + 1);
   hEx0->SetLineWidth(2); hEx1->SetLineWidth(2); hEx2->SetLineWidth(2);

   double sxy = 0, sx = 0, sy = 0, sxx = 0, syy = 0, ssum = 0, ssum2 = 0;
   long nv = 0, nTruth = 0, nReco = 0, nMatch = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      if (tg.IsNull())
         continue;
      TString fs = dir + "/" + tg + "_sim.root", fr = dir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) {
         printf("  skip %-12s (missing)\n", tg.Data());
         continue;
      }
      TH1D *hSel = tg.BeginsWith("gs") ? hEx0 : (tg.BeginsWith("360") ? hEx1 : hEx2);

      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tr = Fr ? (TTree *)Fr->Get("cbmsim") : nullptr;
      if (!ts || !tr) { if (Fs) Fs->Close(); if (Fr) Fr->Close(); continue; }
      TClonesArray *mc = nullptr, *pa = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);

      Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         tr->GetEntry(i);

         // ---- truth: the primary deuteron
         double thTrue = -1, Ttrue = 0, zTrue = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0) break;
            thTrue = std::acos(pz / pp) * TMath::RadToDeg();
            Ttrue = std::sqrt(pp * pp + M_e * M_e) - M_e;
            zTrue = p->GetStartZ(); // cm
            break;
         }
         if (thTrue < 0) continue;
         ++nTruth;
         hT->Fill(thTrue, Ttrue);

         if (!pa || !pa->GetEntriesFast()) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;

         // ---- reconstructed: best angular match, exactly as the gate macros match
         double bd = 1e9, bTh = 0, bT = 0, bZ = 0;
         bool got = false;
         for (auto &track : pe->GetTrackCand()) {
            auto res = spy.Estimate(const_cast<AtTrack &>(track));
            if (!res.valid) continue;
            double th = 180.0 - res.polar * TMath::RadToDeg();
            double p = c * res.brho; // Brho [T*m] -> MeV/c, Z = 1
            double T = std::sqrt(p * p + M_e * M_e) - M_e;
            double d = std::fabs(th - thTrue);
            if (d < bd) { bd = d; bTh = th; bT = T; bZ = res.vertex.Z() / 10.0; got = true; }
         }
         if (!got) continue;
         ++nReco;
         hR->Fill(bTh, bT);
         hV->Fill(zTrue, bZ);
         sx += zTrue; sy += bZ; sxy += zTrue * bZ; sxx += zTrue * zTrue; syy += bZ * bZ; ++nv;
         ssum += zTrue + bZ; ssum2 += (zTrue + bZ) * (zTrue + bZ);
         if (bd > dThetaMax) continue;
         ++nMatch;

         // ---- pre-fit Ex, two-body inversion with the beam energy at the vertex
         // mirror the reconstructed vertex back into the true drift sense (see the header)
         double zUse = driftLength - bZ;
         double Tb = Tb0 - dEdz * zUse;
         if (Tb < 50 || Tb > Tb0 + 20) continue; // vertex outside the chamber -> no sensible beam energy
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double Ed = bT + M_e, pd = std::sqrt(bT * (bT + 2 * M_e));
         double th = bTh * TMath::DegToRad();
         double ER = Eb + M_t - Ed;
         double pRz = pb - pd * std::cos(th), pRt = pd * std::sin(th);
         double m2 = ER * ER - pRz * pRz - pRt * pRt;
         if (m2 <= 0) continue;
         double Ex = std::sqrt(m2) - M_R;
         hEx->Fill(Ex);
         hSel->Fill(Ex);
      }
      printf("  %-12s truth %6ld   reco %6ld   matched %6ld\n", tg.Data(), nTruth, nReco, nMatch);
      Fs->Close();
      Fr->Close();
   }
   delete ta;

   // ---- handedness verdict, printed before it is used anywhere
   double r = (nv * sxy - sx * sy) / std::sqrt((nv * sxx - sx * sx) * (nv * syy - sy * sy));
   double mSum = ssum / std::max(1L, nv);
   double rmsSum = std::sqrt(std::max(0.0, ssum2 / std::max(1L, nv) - mSum * mSum));
   printf("\n  vertex z correlation reco vs truth: r = %+.3f  -> %s\n", r,
          r > 0 ? "SAME sense, no flip applied" : "REVERSED, mirrored back before use");
   printf("  z_true + z_reco = %.2f +- %.2f cm   (drift length used: %.1f cm)%s\n", mSum, rmsSum, driftLength,
          std::fabs(mSum - driftLength) > 2.0 ? "   <-- MISMATCH, the Ex scale below is not trustworthy" : "");

   // ---- analytic loci: solid at mid-chamber, dashed at the entrance and exit energies
   const double TbIn = 598.0, TbOut = 598.0 - 0.957 * driftLength, TbMid = 0.5 * (TbIn + TbOut);
   struct St { double ex; int col; const char *lab; };
   std::vector<St> states = {{0.0, kBlack, "g.s. 1/2^{+}"}, {0.360, kAzure + 2, "0.36 MeV 3/2^{+}"},
                             {2.020, kRed + 1, "2.02 MeV 7/2^{-}"}};
   auto overlay = [&](bool withLegend) {
      auto *lg = new TLegend(0.50, 0.60, 0.88, 0.88);
      lg->SetFillStyle(0);
      lg->SetBorderSize(0);
      for (auto &st : states) {
         auto *gm = ar46_locus(st.ex, TbMid, st.col, 1);
         gm->Draw("L same");
         ar46_locus(st.ex, TbIn, st.col, 2)->Draw("L same");
         ar46_locus(st.ex, TbOut, st.col, 2)->Draw("L same");
         if (withLegend) lg->AddEntry(gm, st.lab, "l");
      }
      if (withLegend) {
         lg->SetHeader("solid: mid-chamber   dashed: entrance / exit");
         lg->Draw();
      }
   };

   TCanvas *c1 = new TCanvas("cK", "pre-fit kinematics", 1400, 1000);
   c1->Divide(2, 2);
   c1->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13); hT->Draw("colz"); overlay(true);
   c1->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.13); hR->Draw("colz"); overlay(false);
   c1->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.13); hV->Draw("colz");
   c1->cd(4);
   hEx0->SetTitle("pre-fit E_{x}(^{47}K), by generated state;E_{x} [MeV];tracks");
   hEx0->Draw("hist"); hEx1->Draw("hist same"); hEx2->Draw("hist same");
   auto *leg = new TLegend(0.62, 0.70, 0.88, 0.88);
   leg->AddEntry(hEx0, "g.s. 1/2^{+}", "l");
   leg->AddEntry(hEx1, "0.36 MeV 3/2^{+}", "l");
   leg->AddEntry(hEx2, "2.02 MeV 7/2^{-}", "l");
   leg->Draw();
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c1->SaveAs(png);

   // FWHM BY GAUSSIAN FIT, not by walking bins from the maximum. The bin-walk estimator this
   // replaced truncates on a single downward fluctuation next to the peak, which made it report
   // 0.765 MeV and 1.800 MeV for two histograms that differ only by a CONSTANT shift of every
   // entry -- a shift cannot change a width, so the estimator was measuring noise. Smoothed
   // half-height crossings set the fit window; the fit itself runs on the unsmoothed histogram.
   auto stat = [&](TH1D *h, const char *lab, double truth) {
      if (h->GetEntries() < 50) { printf("  %-10s too few entries\n", lab); return; }
      TH1D *hs = (TH1D *)h->Clone(TString(h->GetName()) + "_s");
      hs->Smooth(3);
      int b = hs->GetMaximumBin();
      double peak = hs->GetBinCenter(b), half = hs->GetMaximum() / 2, lo = peak, hi = peak;
      for (int k = b; k > 1 && hs->GetBinContent(k) > half; --k) lo = hs->GetBinCenter(k);
      for (int k = b; k < hs->GetNbinsX() && hs->GetBinContent(k) > half; ++k) hi = hs->GetBinCenter(k);
      double w = std::max(0.15, hi - lo);
      TF1 g("g", "gaus", peak - w, peak + w);
      g.SetParameters(h->GetMaximum(), peak, w / 2.355);
      h->Fit(&g, "QNR");
      double mu = g.GetParameter(1), fw = 2.3548 * std::fabs(g.GetParameter(2));
      double qp[3] = {0.25, 0.50, 0.75}, qv[3];
      h->GetQuantiles(3, qv, qp);
      printf("  %-10s median %+6.3f (true %+.2f)   IQR %.3f MeV   gaus-peak %+6.3f FWHM %.3f   entries %.0f\n", lab,
             qv[1], truth, qv[2] - qv[0], mu, fw, h->GetEntries());
      delete hs;
   };
   printf("\n  PRE-FIT Ex, no fitter involved:\n");
   stat(hEx0, "g.s. 1/2+", 0.0);
   stat(hEx1, "0.36 3/2+", 0.36);
   stat(hEx2, "2.02 7/2-", 2.02);
   printf("\n  wrote %s\n\n", png.Data());
}
